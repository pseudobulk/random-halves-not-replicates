"""Model pilot: biolord (0.0.3) on the donor-1 resting / stimulated subset (pilot_common.build_adata; same cells and
split as pilot_cpa.py). Categorical attributes: condition (guide label) and state. Module and trainer settings, n_latent
and the training call are those of biolord's own tutorial (docs/tutorials/biolord_pipeline.ipynb); no tuning.
Counterfactual (as compute_prediction_adata, restricted to the test guides): stimulated control cells from training
(biolord keeps a learned latent per training cell) are given each test guide's condition code and decoded; a guide's
predicted effect is the mean over these cells minus biolord's prediction for the same cells under the control label.
Adapted variant (env BIOLORD_CAT_DIM, BIOLORD_NO_ES=1): the tutorial gives each categorical attribute a 4-dimensional
embedding and stops on a plateau of its biolord metric; with ~2,500 guide categories that run did not fit the guides'
resting effects (fit check r = 0.025), so the variant widens the embedding and trains a fixed number of epochs. Any
variant is reported as adapted, next to the tutorial run.
Usage: CUDA_VISIBLE_DEVICES=<gpu> pilot_biolord.py <data_dir> <out_dir> [max_epochs] [n_source]"""
import os, sys, json, time
import numpy as np, torch
import biolord
from pilot_common import build_adata, effects, score

MODULE = {"decoder_width": 1024, "decoder_depth": 4, "attribute_nn_width": 512, "attribute_nn_depth": 2,
          "n_latent_attribute_categorical": 4, "gene_likelihood": "normal", "reconstruction_penalty": 1e2,
          "unknown_attribute_penalty": 1e1, "unknown_attribute_noise_param": 1e-1, "attribute_dropout_rate": 0.1,
          "use_batch_norm": False, "use_layer_norm": False, "seed": 42}
TRAINER = {"n_epochs_warmup": 0, "latent_lr": 1e-4, "latent_wd": 1e-4, "decoder_lr": 1e-4, "decoder_wd": 1e-4,
           "attribute_nn_lr": 1e-2, "attribute_nn_wd": 4e-8, "step_size_lr": 45, "cosine_scheduler": True,
           "scheduler_final_lr": 1e-5}


def main():
    d, out = sys.argv[1:3]
    maxep = int(sys.argv[3]) if len(sys.argv) > 3 else 500; nsrc = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
    os.makedirs(out, exist_ok=True); torch.set_float32_matmul_precision("high"); t0 = time.time()
    des, adata, genes = build_adata(d)
    adata.obs["condition"] = adata.obs["condition"].astype("category"); adata.obs["state"] = adata.obs["state"].astype("category")
    print("cells", adata.n_obs, "genes", adata.n_vars, dict(adata.obs["split"].value_counts()), "%.0fs" % (time.time() - t0), flush=True)
    biolord.Biolord.setup_anndata(adata, ordered_attributes_keys=None, categorical_attributes_keys=["condition", "state"])
    module = dict(MODULE)
    if os.environ.get("BIOLORD_CAT_DIM"):
        module["n_latent_attribute_categorical"] = int(os.environ["BIOLORD_CAT_DIM"])
    es = os.environ.get("BIOLORD_NO_ES") != "1"
    model = biolord.Biolord(adata=adata, n_latent=32, model_name="cd4_pilot", module_params=module,
                            train_classifiers=False, split_key="split", train_split="train", valid_split="valid", test_split="ood")
    model.train(max_epochs=maxep, batch_size=512, plan_kwargs=TRAINER, early_stopping=es, early_stopping_patience=20,
                check_val_every_n_epoch=10, num_workers=1, enable_checkpointing=False)
    train_s = time.time() - t0
    obs = adata.obs
    src = np.where((obs["condition"] == "ctrl") & (obs["state"] == "stim") & (obs["split"] == "train"))[0]
    src = np.random.default_rng(42).choice(src, min(nsrc, len(src)), replace=False)
    ds = model.get_dataset(adata[src].copy())
    ds = {k: v.to(model.device) for k, v in ds.items()}
    held = (obs["split"] == "ood").to_numpy(); lab = obs["condition"].to_numpy()
    tests = sorted(set(lab[held]))
    # each guide's condition code, read from one of its resting (training) cells
    code = {}
    for g in tests + ["ctrl"]:
        i = np.where((lab == g) & (obs["state"] == "rest").to_numpy())[0][0]
        code[g] = model.get_dataset(adata[[i]].copy())["condition"][0].to(model.device)
    n = len(src)

    @torch.no_grad()
    def predict_as(g):
        dd = dict(ds); dd["condition"] = code[g].unsqueeze(0).repeat(n, *([1] * code[g].ndim))
        mu, _ = model.module.get_expression(dd)
        return mu.mean(0).cpu().numpy()

    model.module.eval()
    base = predict_as("ctrl")
    pred = {g: predict_as(g) - base for g in tests}
    _, res = effects(d)
    sc = score(pred, des, res, genes)
    # fit check: the same guides decoded from resting control cells, against their observed resting effects (training
    # data); a model that cannot reproduce these has not learned per-guide effects at all
    src_r = np.where((obs["condition"] == "ctrl") & (obs["state"] == "rest") & (obs["split"] == "train"))[0]
    src_r = np.random.default_rng(43).choice(src_r, min(nsrc, len(src_r)), replace=False)
    ds = {k: v.to(model.device) for k, v in model.get_dataset(adata[src_r].copy()).items()}; n = len(src_r)
    base_r = predict_as("ctrl")
    R = res["D1_Rest"]; ri = {l: i for i, l in enumerate(R["labels"])}
    fit = [np.corrcoef(predict_as(g) - base_r, R["eff"][ri[g]])[0, 1] for g in tests if g in ri]
    fit_check = dict(resting_fit_pearson=float(np.nanmean(fit)), n=len(fit))
    print("fit check", fit_check, flush=True)
    json.dump(dict(model="biolord", biolord_version=getattr(biolord, "__version__", "?"), train_seconds=train_s,
                   n_source=int(n), scores=sc, fit_check=fit_check, module=module, early_stopping=es, max_epochs=maxep,
                   rules=__doc__),
              open(os.path.join(out, "biolord_scores.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(out, "biolord_preds.npz"), labels=np.array(tests), pred=np.stack([pred[g] for g in tests]))
    print(json.dumps(sc), "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
