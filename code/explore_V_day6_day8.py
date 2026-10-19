"""Explore V: biological-replicate calibration of the cross-line interaction. K562 day 6 (Replogle essential screen)
vs K562 day 8 (Replogle genome-wide screen): same line, lab and platform, different library, time point and batch.
Runs the same gamma pipeline as F26 / F28 / explore_U on this same-line pair, so the "pseudo-gamma" here is what
library + time + batch produce without a cell-line difference. Compare with the lane-split floor (Z1 c: <= 0.3%) and
the real pairs (K562-RPE1 gamma share 0.35, Jurkat-HepG2 0.30).
Self-contained for the A100 storage server (numpy / h5py / scipy only); functions copied unchanged from
make_figs.Line.prep, make_pseudobulk (readers), explore_L_cellcycle, explore_P, explore_R, explore_S, explore_U.
Usage: explore_V_day6_day8.py <day6.h5ad> <day8.h5ad> <outdir>
Inputs: day 6 = wei2025 Replogle_K562essential.h5ad (layers/counts, dense); day 8 = figshare
K562_gwps_raw_singlecell_01.h5ad. Only cells of perturbations with >= 20 cells in both files, plus controls, are
read (chunked). Both files are restricted to their common genes BEFORE CP10K normalisation, so normalisation is on
the same gene set.
Limits: one same-line pair; the two screens also differ in cell number per perturbation (handled by denoising,
which corrects in expectation); exploratory."""
import sys, json, time, numpy as np, h5py

MIN_CELLS = 20
CTRL_NAMES = ("non-targeting", "CTRL", "control")
CHUNK = 50000
S_GENES = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1",
           "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3",
           "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1",
           "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"]
G2M = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO",
       "CENPF", "TACC3", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
       "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
       "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3",
       "CBX5", "CENPA"]
PRO = S_GENES + G2M
P53 = ["MDM2", "FDXR", "TP53I3", "RRM2B", "ZMAT3", "BAX", "CDKN1A", "GDF15", "SESN1", "SESN2", "TNFRSF10B", "AEN",
       "PLK3", "TP53INP1", "DDB2", "XPC", "RPS27L", "TRIAP1", "FAS", "BBC3", "PHLDA3", "SPATA18", "POLH", "TIGAR",
       "CCNG1", "APAF1", "EDA2R", "SULF2"]


# ---- readers (make_pseudobulk) ----
def dec(v):
    return v.decode() if isinstance(v, bytes) else v


def read_codes_cats(f, key):
    obs = f["obs"]; node = obs[key]
    if isinstance(node, h5py.Group):
        codes = node["codes"][:]; cats = np.array([dec(c) for c in node["categories"][:]])
    else:
        codes = node[:]
        if "__categories" in obs and key in obs["__categories"]:
            cats = np.array([dec(c) for c in obs["__categories"][key][:]])
        else:
            cats, codes = np.unique(np.array([dec(v) for v in codes]), return_inverse=True)
    return codes.astype(np.int64), cats


def var_gene_names(f):
    var = f["var"]
    if "gene_name" in var:
        node = var["gene_name"]
        if isinstance(node, h5py.Group):
            cats, codes = node["categories"][:], node["codes"][:]
            return np.array([dec(cats[c]) for c in codes])
        if "__categories" in var and "gene_name" in var["__categories"]:
            cats, codes = var["__categories"]["gene_name"][:], node[:]
            return np.array([dec(cats[c]) for c in codes])
        return np.array([dec(v) for v in node[:]])
    for k in ("_index", "index", "gene_id"):
        if k in var:
            return np.array([dec(v) for v in var[k][:]])
    raise KeyError("no gene names")


def obs_info(path):
    with h5py.File(path) as f:
        codes, cats = read_codes_cats(f, "gene")
        gc, gcats = read_codes_cats(f, "gem_group") if isinstance(f["obs"]["gem_group"], h5py.Group) or \
            "__categories" in f["obs"] and "gem_group" in f["obs"]["__categories"] else (f["obs"]["gem_group"][:], None)
        genes = var_gene_names(f)
    gem = np.asarray(gc if gcats is None else gcats[gc])
    return codes, list(cats), gem, genes


def read_rows(path, rows, cols):
    """Dense rows x cols (float32) from layers/counts or X (dense or CSR), chunked over sorted rows."""
    with h5py.File(path) as f:
        M = f["layers"]["counts"] if "layers" in f and "counts" in f["layers"] else f["X"]
        out = np.zeros((len(rows), len(cols)), np.float32)
        N = M.shape[0] if isinstance(M, h5py.Dataset) else len(M["indptr"]) - 1
        pos = 0
        for s in range(0, N, CHUNK):
            e = min(N, s + CHUNK)
            sel = rows[(rows >= s) & (rows < e)]
            if not len(sel):
                continue
            if isinstance(M, h5py.Dataset):
                blk = M[s:e][:, cols]
            else:
                ip = M["indptr"][s:e + 1]
                data, ind = M["data"][ip[0]:ip[-1]], M["indices"][ip[0]:ip[-1]]
                full = np.zeros((e - s, M.attrs["shape"][1] if "shape" in M.attrs else int(ind.max()) + 1), np.float32)
                rr = np.repeat(np.arange(e - s), np.diff(ip))
                full[rr, ind] = data
                blk = full[:, cols]
            out[pos:pos + len(sel)] = blk[sel - s]; pos += len(sel)
            print("  rows", e, "/", N, flush=True)
    return out


# ---- make_figs.Line.prep ----
class Line:
    def __init__(self, X, codes, cats, genes, gem, ctrl_name):
        self.codes, self.cats = codes, cats
        self.g2j = {g: j for j, g in enumerate(genes)}
        umi = X.sum(1)
        lin = (X / umi[:, None] * 1e4).astype(np.float32)
        self.Y = np.log1p(lin)
        ic = self.cats.index(ctrl_name); isc = self.codes == ic
        ug, gi = np.unique(gem, return_inverse=True)
        gmu = self.Y[isc].mean(0)
        for k in range(len(ug)):
            m = isc & (gi == k)
            if m.sum() < 5:
                m = isc
            self.Y[gi == k] += gmu - self.Y[m].mean(0)
        self.ctrl = np.where(isc)[0]

    def n_cells(self, t):
        return int((self.codes == self.cats.index(t)).sum())


# ---- explore_L_cellcycle ----
def split(genes, gpos, rng):
    g = [x for x in genes if x in gpos]; g = list(rng.permutation(g)); return g[: len(g) // 2], g[len(g) // 2:]


def phases(Y, mu_ctrl, set_ix, rng):
    bins = np.digitize(mu_ctrl, np.quantile(mu_ctrl, np.linspace(0, 1, 26)[1:-1]))
    ctl = np.concatenate([rng.choice(np.where(bins == bins[j])[0], 10, replace=False) for j in set_ix])
    return Y[:, set_ix].mean(1) - Y[:, ctl].mean(1)


def line_parts(L, perts, gi, s1, g1, seed):
    rng = np.random.default_rng(seed)
    Y = L.Y[:, gi]; mu = Y[L.ctrl].mean(0)
    sS, sG = phases(Y, mu, s1, rng), phases(Y, mu, g1, rng)
    ph = np.where((sS <= 0) & (sG <= 0), 0, np.where(sS >= sG, 1, 2))
    c = L.ctrl
    fc = np.array([(ph[c] == k).mean() for k in range(3)])
    mph = np.stack([Y[c][ph[c] == k].mean(0) for k in range(3)])
    P, G = len(perts), len(gi)
    act, pred = np.zeros((2, P, G), np.float32), np.zeros((2, P, G), np.float32)
    r2 = np.random.default_rng(seed + 1)
    for k, p in enumerate(perts):
        rows = r2.permutation(np.where(L.codes == L.cats.index(p))[0]); h = len(rows) // 2
        for hh, rr in enumerate((rows[:h], rows[h:2 * h])):
            f = np.array([(ph[rr] == q).mean() for q in range(3)])
            act[hh, k] = Y[rr].mean(0) - mu
            pred[hh, k] = (f - fc) @ mph
    return dict(fc=fc, act=act, pred=pred)


def tmpl(XA, XB):
    X = np.stack([XA, XB]); m = X.mean((0, 1)); return X - m - (X.mean(1) - m)[:, None, :]


# ---- explore_P / R / S ----
def share(g, gA, F):
    sig = (g ** 2).sum() - ((gA - g) ** 2).sum() / 2
    return float(((g ** 2).sum() - ((g - F) ** 2).sum()) / sig)


def gain_fit(res, beta, cols):
    F = np.zeros_like(res); b = beta[:, cols]
    k = (res[:, cols] * b).sum(0) / np.maximum((b ** 2).sum(0), 1e-12)
    F[:, cols] = b * k
    return F, k


def wshare(W, RA, RB):
    cov = W.sum(1) > 0
    Wn = W[cov] / W[cov].sum(1, keepdims=True)
    P_ = Wn @ RA
    T, D = RB[cov], (RA - RB)[cov]
    sig = (T ** 2).sum() - (D ** 2).sum() / 2
    c = (T * P_).sum() / max((P_ ** 2).sum(), 1e-12)
    return dict(share_covered=float(((T ** 2).sum() - ((T - c * P_) ** 2).sum()) / sig))


def cosine_knn(XA, k):
    n = XA / np.maximum(np.linalg.norm(XA, axis=1, keepdims=True), 1e-12)
    S = n @ n.T; np.fill_diagonal(S, -np.inf)
    W = np.zeros_like(S); np.put_along_axis(W, np.argsort(-S, 1)[:, :k], 1.0, 1)
    return W


# ---- explore_U.pair_stats ----
def pair_stats(LA, LB):
    genes = sorted(set(LA.g2j) & set(LB.g2j)); gpos = {g: k for k, g in enumerate(genes)}
    rng = np.random.default_rng(7)
    s1, _ = split(S_GENES, gpos, rng); g1, _ = split(G2M, gpos, rng)
    i_s1, i_g1 = np.array([gpos[g] for g in s1]), np.array([gpos[g] for g in g1])
    excl = set(P53) | set(PRO) | {"TP53"}
    okp = [{c for c in L.cats if c not in CTRL_NAMES and L.n_cells(c) >= MIN_CELLS} for L in (LA, LB)]
    perts = sorted((okp[0] & okp[1]) - excl)
    giA, giB = np.array([LA.g2j[g] for g in genes]), np.array([LB.g2j[g] for g in genes])
    mA = LA.Y[LA.ctrl][:, giA].mean(0); mB = LB.Y[LB.ctrl][:, giB].mean(0)
    A = line_parts(LA, perts, giA, i_s1, i_g1, 1)
    Bp = line_parts(LB, perts, giB, i_s1, i_g1, 2)
    for arr in (A["act"], A["pred"], Bp["act"], Bp["pred"]):
        for k, t in enumerate(perts):
            if t in gpos:
                arr[:, k, gpos[t]] = 0
    phase_call = set(s1) | set(g1)
    ev = np.array([gpos[x] for x in genes if x not in phase_call])
    f64 = lambda a: a[:, ev].astype(np.float64)
    TA, TB = tmpl(A["act"][0], Bp["act"][0]), tmpl(A["act"][1], Bp["act"][1])
    gA, g = f64(TA[0] - TA.mean(0)), f64(TB[0] - TB.mean(0))
    bA, bB = f64(TA.mean(0)), f64(TB.mean(0))
    pc = tmpl(A["pred"][0], Bp["pred"][0]); comp = f64(pc[0] - pc.mean(0))
    base = bA * ((mA - mB) / 2)[ev]
    out = dict(n_perturbations=len(perts), n_genes=int(len(ev)), phase_frac_ctrl=[A["fc"].tolist(), Bp["fc"].tolist()])
    bs, gs = float((bA * bB).sum()), float((gA * g).sum())
    out["beta_signal"], out["gamma_signal"] = bs, gs
    out["gamma_share_of_specific"] = gs / (bs + gs)
    out["gamma_signal_frac_of_raw"] = float(((g ** 2).sum() - ((gA - g) ** 2).sum() / 2) / (g ** 2).sum())
    c1 = np.linalg.lstsq(comp.ravel()[:, None], g.ravel(), rcond=None)[0]
    X2 = np.column_stack([comp.ravel(), base.ravel()])
    F2 = (X2 @ np.linalg.lstsq(X2, g.ravel(), rcond=None)[0]).reshape(g.shape)
    F4 = F2 + gain_fit(g - F2, bA, np.arange(len(ev)))[0]
    out["nested"] = dict(comp=share(g, gA, comp * c1[0]), comp_base=share(g, gA, F2), comp_base_gain=share(g, gA, F4))
    out["borrow_raw_gamma_k10"] = dict(ideal=wshare(cosine_knn(gA, 10), gA, g)["share_covered"],
                                       picked_by_beta=wshare(cosine_knn(bA, 10), gA, g)["share_covered"])
    r2 = np.random.default_rng(5)
    pa, pb = np.empty_like(gA), np.empty_like(g)
    for i in range(len(gA)):
        o = r2.permutation(gA.shape[1]); pa[i], pb[i] = gA[i, o], g[i, o]
    out["borrow_raw_gamma_k10"]["floor_genes_shuffled"] = wshare(cosine_knn(pa, 10), pa, pb)["share_covered"]
    # where does the same-line difference live? gamma energy on proliferation / p53 / ribosomal / mito genes
    gev = [genes[j] for j in ev]
    groups = dict(proliferation=set(PRO), p53=set(P53),
                  ribosomal={x for x in gev if x.startswith(("RPL", "RPS"))}, mito={x for x in gev if x.startswith("MT-")})
    gsig = (gA * g).sum(0)
    out["gamma_signal_by_gene_group"] = {nm: float(gsig[[j for j, x in enumerate(gev) if x in s]].sum() / gs)
                                         for nm, s in groups.items()}
    return out


def load(path, keep_names):
    codes, cats, gem, genes = obs_info(path)
    ctrl = [c for c in CTRL_NAMES if c in cats][0]
    keep_codes = {cats.index(c) for c in keep_names | {ctrl} if c in cats}
    rows = np.where(np.isin(codes, list(keep_codes)))[0]
    return dict(codes=codes, cats=cats, gem=gem, genes=genes, ctrl=ctrl, rows=rows)


def main():
    p6, p8, out = sys.argv[1], sys.argv[2], sys.argv[3]
    t0 = time.time()
    meta = {}
    for nm, p in (("d6", p6), ("d8", p8)):
        codes, cats, gem, genes = obs_info(p)
        n = np.bincount(codes, minlength=len(cats))
        meta[nm] = dict(cats=cats, n=n, genes=genes)
        print(nm, "cells", len(codes), "perts", len(cats), "genes", len(genes), flush=True)
    ok = [{c for c, k in zip(m["cats"], m["n"]) if k >= MIN_CELLS} for m in meta.values()]
    keep = (ok[0] & ok[1]) - set(CTRL_NAMES)
    common = sorted(set(meta["d6"]["genes"]) & set(meta["d8"]["genes"]))
    print("shared perts", len(keep), "common genes", len(common), flush=True)
    lines = {}
    for nm, p in (("d6", p6), ("d8", p8)):
        d = load(p, keep)
        gidx = {g: j for j, g in enumerate(d["genes"])}
        cols = np.array([gidx[g] for g in common])
        X = read_rows(p, d["rows"], cols)
        sub_codes = d["codes"][d["rows"]]
        lines[nm] = Line(X, sub_codes, d["cats"], common, d["gem"][d["rows"]], d["ctrl"])
        print(nm, "loaded", X.shape, "ctrl", len(lines[nm].ctrl), round(time.time() - t0), "s", flush=True)
        del X
    res = {"rules": __doc__, "n_shared_perts_read": len(keep), "n_common_genes": len(common)}
    res["day6_vs_day8"] = d = pair_stats(lines["d6"], lines["d8"])
    print("day6_vs_day8", json.dumps(d, default=lambda v: round(float(v), 4) if np.isscalar(v) else v), flush=True)
    json.dump(res, open(out + "/metrics_explore_V_day6_day8.json", "w"), indent=1,
              default=lambda v: float(v) if np.isscalar(v) else str(v))
    print("saved", round(time.time() - t0), "s", flush=True)


if __name__ == "__main__":
    main()
