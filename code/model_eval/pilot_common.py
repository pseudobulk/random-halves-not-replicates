"""Model pilot, shared code (notes/model_pilot_shortlist.md). Loads the donor-1 subsets written by
code/a100/cd4_pilot_extract.py, builds guide-level effects and scores predictions of stimulated effects against two
truths for each test guide g of target t:
  A  same clones          the stimulated cells of guide g itself (its resting cells, the input, share its clones)
  B  independent replicate the stimulated cells of the other guide of t (separate transduction, no shared clone)
The input and truth cells are different cells, so sampling noise does not link them; A minus B is what a prediction
gains from the clones (plus, for expressed targets, any guide-specific biology; the strength-matched set limits that).
Expression: log1p CP10K with the cell's UMI over all genes. Effect: mean over a guide's cells minus mean over control
cells of the same state. Labels: '<target>|g<gp>'; controls 'ctrl'. Donor: env PILOT_DONOR (default D1; model_eval
step 1 runs D1-D4; env PILOT_STATES=K562,RPE1 for the negative control)."""
import os, json
import numpy as np
from scipy.sparse import csr_matrix

DONOR = os.environ.get("PILOT_DONOR", "D1")
STATES = (DONOR + "_Rest", DONOR + "_Stim48hr")
if os.environ.get("PILOT_STATES"):   # source,target files of another pair, e.g. K562,RPE1 (model_eval negative control)
    DONOR = os.environ["PILOT_STATES"].replace(",", "_to_")
    STATES = tuple(os.environ["PILOT_STATES"].split(","))


def load(d, st):
    z = np.load(os.path.join(d, "pilot_%s.npz" % st), allow_pickle=False)
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    lab = np.array(["ctrl" if k == "control" else "%s|g%d" % (t, p) for t, p, k in zip(z["target"], z["gp"], z["kind"])])
    return dict(X=X, umi=z["umi"], lab=lab, kind=z["kind"], target=z["target"], genes=[str(g) for g in z["genes"]])


def lognorm(X, umi):
    L = X.astype(np.float32).tocsr(copy=True)
    L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
    return L


def group_means(L, lab, labels):
    """dense mean expression per label (rows of L that carry it)"""
    idx = {l: i for i, l in enumerate(labels)}
    rows = np.array([idx.get(l, -1) for l in lab]); ok = rows >= 0
    M = csr_matrix((np.ones(ok.sum(), np.float32), (rows[ok], np.where(ok)[0])), shape=(len(labels), L.shape[0]))
    n = np.asarray(M.sum(1)).ravel()
    out = np.asarray((M @ L).todense()) / np.maximum(n, 1)[:, None]
    return out.astype(np.float32), n


def design(d):
    return json.load(open(os.path.join(d, "pilot_design.json")))


def kinds(des):
    """test groups: inert and strong, plus any listed in extra_kinds (model_eval negative control)"""
    return ("inert", "strong") + tuple(des.get("extra_kinds", []))


def test_targets(des):
    return {t for k in kinds(des) for t in des[k]}


def build_adata(d):
    """one AnnData of both states for the cell-level models: log1p CP10K, obs condition (guide label / 'ctrl'),
    state ('rest' / 'stim'), split: stimulated cells of test guides 'ood' (held out); all other cells 'train', 5% of
    them at random (seed 0) 'valid'. The resting cells of test guides stay in training (context transfer)."""
    import pandas as pd, anndata as ad
    des = design(d); test = test_targets(des)
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


def effects(d):
    """guide-level effects per state; returns dict with labels, R (rest effects), S (stim effects), control means"""
    des = design(d); res = {}
    for st in STATES:
        D = load(d, st); L = lognorm(D["X"], D["umi"])
        labels = sorted(set(D["lab"]) - {"ctrl"})
        M, n = group_means(L, D["lab"], labels + ["ctrl"])
        res[st] = dict(labels=labels, eff=M[:-1] - M[-1], ctrl=M[-1], n=n[:-1], genes=D["genes"])
    assert res[STATES[0]]["genes"] == res[STATES[1]]["genes"]
    return des, res


def test_pairs(des):
    """(target, kind, guide label, other guide label) for every test guide"""
    out = []
    for kind in kinds(des):
        for t in des[kind]:
            for g in (0, 1):
                out.append((t, kind, "%s|g%d" % (t, g), "%s|g%d" % (t, 1 - g)))
    return out


def score(pred, des, res, genes, gene_mask=None):
    """pred: dict guide label -> predicted stimulated effect (vector over genes). Pearson over genes, the target gene
    excluded (and genes outside gene_mask, for models whose gene axis misses some pilot genes); mean over test guides
    per kind, for truth A (same guide) and B (other guide)."""
    S = res[STATES[1]]; si = {l: i for i, l in enumerate(S["labels"])}; gi = {g: j for j, g in enumerate(genes)}
    rows = {k: {"A": [], "B": []} for k in kinds(des)}
    for t, kind, g, o in test_pairs(des):
        if g not in pred or g not in si or o not in si:
            continue
        keep = np.ones(len(genes), bool) if gene_mask is None else np.asarray(gene_mask, bool).copy()
        if t in gi:
            keep[gi[t]] = False
        p = pred[g][keep]
        for k, lab in (("A", g), ("B", o)):
            y = S["eff"][si[lab]][keep]
            rows[kind][k].append(np.corrcoef(p, y)[0, 1] if p.std() > 0 else np.nan)
    out = {}
    for kind, v in rows.items():
        a, b = np.array(v["A"]), np.array(v["B"])
        out[kind] = dict(n=int(len(a)), pearson_A=float(np.nanmean(a)), pearson_B=float(np.nanmean(b)),
                         leak=float(np.nanmean(a - b)), leak_se=float(np.nanstd(a - b) / np.sqrt(max(np.isfinite(a - b).sum(), 1))))
    return out
