"""Model pilot: CPA (cpa-tools 0.8.8) on the donor-1 resting / stimulated subset (pilot_common.py).
Setup: perturbation = guide label, categorical covariate = state, log1p CP10K input (is_count_data False, Gaussian
reconstruction), hyperparameters of scPerturBench's myCPA.py (the published settings; no tuning). No perturbation
embeddings: CPA learns one per guide.
Splits: stimulated cells of test guides = 'ood' (held out); all other cells train, 5% of them at random (seed 0) 'valid'
for early stopping. The resting cells of test guides are in training, as in a context-transfer task.
Counterfactual: the held-out rows are replaced by stimulated control cells drawn at random (seed 42, with replacement)
and predicted with their own guide label and state; a guide's predicted effect is the mean prediction over its rows
minus CPA's prediction for stimulated control cells under the control label. Scored by pilot_common.score.
Usage: CUDA_VISIBLE_DEVICES=<gpu> pilot_cpa.py <data_dir> <out_dir> [max_epochs]"""
import os, sys, json, time
import numpy as np, pandas as pd, anndata as ad, torch
import cpa
import cpa._data as CPA_DATA
from pilot_common import STATES, load, lognorm, design, effects, score

NW = int(os.environ.get("CPA_WORKERS", "0"))
_splitter_init = CPA_DATA.AnnDataSplitter.__init__


def _splitter_init_workers(self, *a, **k):
    # CPA builds its loaders without num_workers (single-core batching); optional workers for speed only
    if NW > 0:
        k.setdefault("num_workers", NW)
        k.setdefault("persistent_workers", True)
    _splitter_init(self, *a, **k)


CPA_DATA.AnnDataSplitter.__init__ = _splitter_init_workers

# CPA computes a logged R2 for every perturbation in every training batch, with a GPU sync per perturbation; with ~2,500
# guide labels this makes a step CPU-bound (~10 min per epoch). The value is only logged during training and does not
# enter the loss, so it is skipped while the module is in training mode; validation (module.eval(), which early
# stopping reads) computes it as usual. Speed only.
import cpa._module as CPA_MODULE
_r2 = CPA_MODULE.CPAModule.r2_metric


def _r2_eval_only(self, *a, **k):
    return (0.0, 0.0) if self.training else _r2(self, *a, **k)


CPA_MODULE.CPAModule.r2_metric = _r2_eval_only

AE = {"n_latent": 128, "recon_loss": "gauss", "doser_type": "linear", "n_hidden_encoder": 512, "n_layers_encoder": 2,
      "n_hidden_decoder": 256, "n_layers_decoder": 2, "use_batch_norm_encoder": True, "use_layer_norm_encoder": False,
      "use_batch_norm_decoder": False, "use_layer_norm_decoder": False, "dropout_rate_encoder": 0.0,
      "dropout_rate_decoder": 0.2, "variational": False, "seed": 1117}
TR = {"n_epochs_kl_warmup": None, "n_epochs_pretrain_ae": 30, "n_epochs_adv_warmup": 100, "n_epochs_mixup_warmup": 10,
      "mixup_alpha": 0.2, "adv_steps": 2, "n_hidden_adv": 128, "n_layers_adv": 2, "use_batch_norm_adv": True,
      "use_layer_norm_adv": False, "dropout_rate_adv": 0.2, "reg_adv": 50.0, "pen_adv": 1.0, "lr": 0.0003, "wd": 4e-07,
      "adv_lr": 0.0003, "adv_wd": 4e-07, "adv_loss": "cce", "doser_lr": 0.0003, "doser_wd": 4e-07,
      "do_clip_grad": False, "gradient_clip_value": 5.0, "step_size_lr": 45}


def build(d):
    des = design(d); test = set(des["inert"]) | set(des["strong"])
    Xs, obs = [], []
    for st, short in zip(STATES, ("rest", "stim")):
        D = load(d, st); L = lognorm(D["X"], D["umi"])
        Xs.append(L.toarray())
        held = np.array([short == "stim" and t in test for t in D["target"]]) & (D["kind"] != "control")
        obs.append(pd.DataFrame({"condition": D["lab"], "state": short, "held": held}))
        genes = D["genes"]
    X = np.vstack(Xs).astype(np.float32); o = pd.concat(obs, ignore_index=True)
    rng = np.random.default_rng(0)
    o["split"] = np.where(o["held"], "ood", np.where(rng.random(len(o)) < 0.05, "valid", "train"))
    o.index = [str(i) for i in range(len(o))]
    return des, ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=genes)), genes


def main():
    d, out = sys.argv[1:3]; maxep = int(sys.argv[3]) if len(sys.argv) > 3 else 500
    os.makedirs(out, exist_ok=True); torch.set_float32_matmul_precision("high"); t0 = time.time()
    des, adata, genes = build(d)
    print("cells", adata.n_obs, "genes", adata.n_vars, dict(adata.obs["split"].value_counts()), "%.0fs" % (time.time() - t0), flush=True)
    cpa.CPA.setup_anndata(adata, perturbation_key="condition", control_group="ctrl", dosage_key=None,
                          categorical_covariate_keys=["state"], is_count_data=False, max_comb_len=1)
    model = cpa.CPA(adata=adata, split_key="split", train_split="train", valid_split="valid", test_split="ood", **AE)
    # adapted variant (env CPA_NO_ES=1): the default run stopped at epoch 40, inside the adversarial warm-up, and did not
    # fit the guides' resting effects (fit check r = 0.003); the variant trains a fixed max_epochs
    patience = 10 ** 6 if os.environ.get("CPA_NO_ES") == "1" else 5
    model.train(max_epochs=maxep, use_gpu=0, batch_size=1024, plan_kwargs=TR, early_stopping_patience=patience,
                check_val_every_n_epoch=5, save_path=os.path.join(out, "model"))
    train_s = time.time() - t0
    held = (adata.obs["split"] == "ood").to_numpy()
    sctrl = ((adata.obs["condition"] == "ctrl") & (adata.obs["state"] == "stim")).to_numpy()
    rng = np.random.default_rng(42)
    a = adata.copy()
    Xn = a.X.copy(); Xn[held] = adata.X[rng.choice(np.where(sctrl)[0], held.sum(), replace=True)]; a.X = Xn
    model.predict(a, batch_size=2048)
    P = np.asarray(a.obsm["CPA_pred"])
    base = P[sctrl].mean(0)   # CPA's own prediction for stimulated controls
    lab = a.obs["condition"].to_numpy()
    pred = {g: P[held & (lab == g)].mean(0) - base for g in np.unique(lab[held])}
    _, res = effects(d)
    sc = score(pred, des, res, genes)
    json.dump(dict(model="CPA", cpa_version=cpa.__version__, train_seconds=train_s, scores=sc, rules=__doc__),
              open(os.path.join(out, "cpa_scores.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(out, "cpa_preds.npz"), labels=np.array(list(pred)), pred=np.stack(list(pred.values())))
    print(json.dumps(sc), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
