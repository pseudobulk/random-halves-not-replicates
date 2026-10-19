"""CD4 T single cells (Marson 2025 genome-scale CRISPRi, GSE314342 / S3 marson2025_data), one donor-condition file per
run, on the A100 storage server (numpy / h5py / scipy only, no anndata). One streaming pass over the CSR matrix; writes
compact aggregates for cd4_sc_analyze.py. Usage: cd4_sc_extract.py <h5ad> <outdir> [workers] [sample_chunks]

Cells: low_quality cells and cells without a perturbed_gene_name are dropped; controls = guide_type 'non-targeting'
(own group, index n_cats). Expression: log1p CP10K (row sums over all genes), as make_figs.Line, without gem centring.
Phase calls as explore_L_cellcycle: Tirosh S / G2M genes split at random (seed 7) into half 1 (phase calls) and half 2;
score = mean log1p CP10K of the half-1 set minus mean of a control set (10 random genes from the same of 25 expression
bins per set gene; bins from the control mean in a random sample of chunks); G1 if both scores <= 0, else the larger.
Per-cell output: group, phase, S / G2M scores, UMI, raw count of the cell's own target gene, raw UMI on all Tirosh S and
on all G2M genes, lane, random half (seed 11).
Aggregates (float32 sums of log1p CP10K over all genes): per (group, phase) and per (group, half), with cell counts.
Control-cell raw-count moments per lane for the escape model of explore_H5 (rate = sum y / sum u per lane and gene;
NB dispersion from sum y^2, sum u*y, sum u^2)."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time, numpy as np, h5py
from multiprocessing import Pool
from scipy.sparse import csr_matrix

S_GENES = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1",
           "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3",
           "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1",
           "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"]
G2M = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO",
       "CENPF", "TACC3", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
       "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
       "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3",
       "CBX5", "CENPA"]
CHUNK = 4000
G_ = {}   # per-process globals (set before fork)


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def cat(o, k):
    return o[k]["codes"][:], dec(o[k]["categories"][:])


def split(genes, gpos, rng):
    g = list(rng.permutation([x for x in genes if x in gpos]))
    return g[: len(g) // 2], g[len(g) // 2:]


def read_rows(f, r0, r1):
    ip = f["X"]["indptr"][r0:r1 + 1]
    a, b = int(ip[0]), int(ip[-1])
    X = csr_matrix((f["X"]["data"][a:b].astype(np.float32), f["X"]["indices"][a:b].astype(np.int32), ip - a),
                   shape=(r1 - r0, G_["G"]))
    X.sum_duplicates()
    return X


def lognorm(X, umi):
    L = X.copy()
    L.data = np.log1p(L.data / np.repeat(umi, np.diff(L.indptr)) * 1e4)
    return L


def sample_ctrl_mean(fn, n_chunks):
    rng = np.random.default_rng(3)
    starts = np.sort(rng.choice(np.arange(0, G_["N"] - CHUNK, CHUNK), n_chunks, replace=False))
    s, n = np.zeros(G_["G"]), 0
    with h5py.File(fn, "r") as f:
        for r0 in starts:
            m = G_["ctrl"][r0:r0 + CHUNK]
            if not m.any():
                continue
            X = read_rows(f, r0, r0 + CHUNK)[np.where(m)[0]]
            umi = np.asarray(X.sum(1)).ravel(); ok = umi > 0
            s += np.asarray(lognorm(X[ok], umi[ok]).sum(0)).ravel(); n += int(ok.sum())
    return s / n, n


def worker(args):
    wid, r0w, r1w = args
    g = G_; Gn, NG = g["G"], g["n_grp"]
    acc_ph = np.zeros((NG * 3, Gn), np.float32); acc_h = np.zeros((NG * 2, Gn), np.float32)
    nl = g["n_lane"]
    cy, cy2, cuy = (np.zeros((nl, Gn)) for _ in range(3)); cu, cu2 = np.zeros(nl), np.zeros(nl)
    n = r1w - r0w
    cell = dict(phase=np.full(n, -1, np.int8), sS=np.zeros(n, np.float32), sG=np.zeros(n, np.float32),
                umi=np.zeros(n, np.float32), y_t=np.zeros(n, np.float32), rawS=np.zeros(n, np.float32),
                rawG=np.zeros(n, np.float32))
    with h5py.File(g["fn"], "r") as f:
        for r0 in range(r0w, r1w, CHUNK):
            r1 = min(r0 + CHUNK, r1w); o = r0 - r0w
            X = read_rows(f, r0, r1)
            umi = np.asarray(X.sum(1)).ravel()
            grp = g["grp"][r0:r1].copy(); grp[umi <= 0] = -1
            v = np.where(grp >= 0)[0]
            if not len(v):
                continue
            X, umi, grp = X[v], umi[v], grp[v]
            sl = o + v
            small = np.log1p(X[:, g["pcols"]].toarray() / umi[:, None] * 1e4)
            sS = small[:, g["iS"]].mean(1) - small[:, g["iSc"]].mean(1)
            sG = small[:, g["iG"]].mean(1) - small[:, g["iGc"]].mean(1)
            ph = np.where((sS <= 0) & (sG <= 0), 0, np.where(sS >= sG, 1, 2)).astype(np.int8)
            rows = np.repeat(np.arange(len(v)), np.diff(X.indptr))
            tc = g["tcol"][r0:r1][v]
            hit = X.indices == tc[rows]
            cell["y_t"][sl] = np.bincount(rows[hit], weights=X.data[hit], minlength=len(v))
            cell["rawS"][sl] = np.asarray(X[:, g["colS"]].sum(1)).ravel()
            cell["rawG"][sl] = np.asarray(X[:, g["colG"]].sum(1)).ravel()
            cell["phase"][sl], cell["sS"][sl], cell["sG"][sl], cell["umi"][sl] = ph, sS, sG, umi
            L = lognorm(X, umi)
            half = g["half"][r0:r1][v]
            for acc, idx in ((acc_ph, grp * 3 + ph), (acc_h, grp * 2 + half)):
                R = (csr_matrix((np.ones(len(v), np.float32), (idx, np.arange(len(v)))), shape=(acc.shape[0], len(v)))
                     @ L).tocoo()
                acc[R.row, R.col] += R.data
            c = np.where(grp == NG - 1)[0]
            if len(c):
                ln = g["lane"][r0:r1][v][c]; Xc = X[c]; uc = umi[c]
                for k in np.unique(ln):
                    m = ln == k; Xk = Xc[np.where(m)[0]]; uk = uc[m]
                    cy[k] += np.asarray(Xk.sum(0)).ravel()
                    cy2[k] += np.asarray(Xk.multiply(Xk).sum(0)).ravel()
                    cuy[k] += np.asarray(Xk.T @ uk).ravel()
                    cu[k] += uk.sum(); cu2[k] += (uk ** 2).sum()
    np.save("%s/part_ph_%03d.npy" % (g["out"], wid), acc_ph)
    np.save("%s/part_h_%03d.npy" % (g["out"], wid), acc_h)
    return wid, r0w, cell, dict(cy=cy, cy2=cy2, cuy=cuy, cu=cu, cu2=cu2)


def main():
    fn, out = sys.argv[1], sys.argv[2]
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    nsamp = int(sys.argv[4]) if len(sys.argv) > 4 else 150
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    with h5py.File(fn, "r") as f:
        o = f["obs"]
        pc, pcats = cat(o, "perturbed_gene_name"); gc, gcats = cat(o, "guide_type"); lc, lcats = cat(o, "lane_id")
        lq = o["low_quality"][:]
        genes = dec(f["var"]["gene_name"][:])
        N, Gn = f["X"].attrs["shape"]
    gpos = {x: k for k, x in enumerate(genes)}
    ctrl_codes = [k for k, x in enumerate(gcats) if "non" in x.lower()]
    ctrl = np.isin(gc, ctrl_codes)
    NG = len(pcats) + 1
    grp = np.where(ctrl, NG - 1, pc).astype(np.int64)
    grp[(pc < 0) & ~ctrl] = -1; grp[lq.astype(bool)] = -1
    t_of_cat = np.array([gpos.get(x, -1) for x in pcats] + [-1])
    tcol = np.where(grp >= 0, t_of_cat[np.maximum(grp, 0)], -1)
    G_.update(fn=fn, out=out, N=int(N), G=int(Gn), ctrl=ctrl & (grp >= 0), grp=grp, tcol=tcol,
              lane=lc.astype(np.int64), n_lane=len(lcats), n_grp=NG,
              half=(np.random.default_rng(11).random(int(N)) < 0.5).astype(np.int64))
    mu0, n_s = sample_ctrl_mean(fn, nsamp)
    print("ctrl mean from %d sampled control cells, %.0fs" % (n_s, time.time() - t0), flush=True)
    rng = np.random.default_rng(7)
    s1, s2 = split(S_GENES, gpos, rng); g1, g2 = split(G2M, gpos, rng)
    bins = np.digitize(mu0, np.quantile(mu0, np.linspace(0, 1, 26)[1:-1]))
    pick = lambda st: np.concatenate([rng.choice(np.where(bins == bins[gpos[x]])[0], 10, replace=False) for x in st])
    sets = dict(S=np.array([gpos[x] for x in s1]), Sc=pick(s1), G=np.array([gpos[x] for x in g1]), Gc=pick(g1))
    pcols = np.unique(np.concatenate(list(sets.values()))); cpos = {c: k for k, c in enumerate(pcols)}
    G_.update(pcols=pcols, **{"i" + k: np.array([cpos[c] for c in v]) for k, v in sets.items()},
              colS=np.array([gpos[x] for x in S_GENES if x in gpos]), colG=np.array([gpos[x] for x in G2M if x in gpos]))
    edges = np.linspace(0, int(N), nw + 1).astype(int)
    jobs = [(w, edges[w], edges[w + 1]) for w in range(nw)]
    cell = {k: None for k in ("phase", "sS", "sG", "umi", "y_t", "rawS", "rawG")}
    for k in cell:
        cell[k] = np.zeros(int(N), np.int8 if k == "phase" else np.float32)
    cell["phase"][:] = -1
    mom = None
    with Pool(nw) as pool:
        for wid, r0, c, m in pool.imap_unordered(worker, jobs):
            for k in cell:
                cell[k][r0:r0 + len(c[k])] = c[k]
            mom = m if mom is None else {k: mom[k] + m[k] for k in m}
            print("worker %d done %.0fs" % (wid, time.time() - t0), flush=True)
    acc_ph = sum(np.load("%s/part_ph_%03d.npy" % (out, w)) for w in range(nw))
    acc_h = sum(np.load("%s/part_h_%03d.npy" % (out, w)) for w in range(nw))
    ok = cell["phase"] >= 0
    n_ph = np.bincount(grp[ok] * 3 + cell["phase"][ok], minlength=NG * 3).reshape(NG, 3)
    n_h = np.bincount(grp[ok] * 2 + G_["half"][ok], minlength=NG * 2).reshape(NG, 2)
    np.save(out + "/sum_phase.npy", acc_ph.reshape(NG, 3, Gn)); np.save(out + "/sum_half.npy", acc_h.reshape(NG, 2, Gn))
    np.savez(out + "/cells.npz", grp=grp, lane=G_["lane"], half=G_["half"], tcol=tcol, **cell)
    np.savez(out + "/ctrl_moments.npz", n_phase=n_ph, n_half=n_h, mu0=mu0, **mom)
    json.dump(dict(file=fn, n_cells=int(N), n_genes=int(Gn), genes=genes, groups=pcats + ["__ctrl__"],
                   lanes=lcats, phase_genes=dict(S1=s1, S2=s2, G1=g1, G2=g2), n_ctrl_sampled=n_s,
                   n_kept=int(ok.sum()), n_ctrl_kept=int(n_ph[-1].sum()), rules=__doc__),
              open(out + "/meta.json", "w"))
    for w in range(nw):
        os.remove("%s/part_ph_%03d.npy" % (out, w)); os.remove("%s/part_h_%03d.npy" % (out, w))
    print("saved %s, kept %d cells, %.0fs" % (out, ok.sum(), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
