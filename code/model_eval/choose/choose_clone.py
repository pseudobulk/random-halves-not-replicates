"""CHOOSE (Li, Fleck et al. 2023, Nature; Zenodo 7083558 object = batch 1, 8 libraries, 28,344 cells with one gRNA
cassette each; clone barcodes from E-MTAB-13148 parsed by cbc_parse.py): the direct test of the clone mechanism and of
clone strength = ICC / K_eff, with clone labels measured per cell.

Clone calls. All (library, cell barcode, UMI, UCB) read counts of a library are pooled over its amplicon runs; only
cell barcodes of the object are kept. As in the paper's Methods: within a cell barcode-UMI pair, UCBs with < 20% of the
pair's top read count are dropped; triples below 1% of the library's median (over cells) top-UMI read count are dropped.
UCBs one substitution away from a more abundant UCB (UMIs summed over all cells) are merged into it. A cell gets the
clone of its top UCB if that UCB has >= MINU UMIs and >= 2x the UMIs of the second. Check: a clone should carry one gRNA.

Expression. log1p CP10K; genes with mean >= 0.1 over all cells. Each cell's deviation = expression minus the mean of
the other cells of its reference group in its library, leaving out the cells of its own perturbation (as for Jin 2020):
  type    reference group = fine cell type (celltype_cl_coarse2): clone state within cell types (main)
  library reference group = the whole library: also counts differences in cell-type composition
Per-cell variance V = summed over genes, of cells around their reference group mean.

Per perturbation (gRNA, all cell types pooled; cells with a clone call; >= MINC cells):
  random   random halves of the cells (NSPLIT splits), cross-half signal sum_g x_A,g x_B,g
  clone    clones assigned at random to two halves (whole clones; NSPLIT splits; >= 10 cells per half)
  library  libraries split into two random groups (separate organoid pools; NSPLIT splits)
  clone part = random - clone; strength = clone part / V (same scale as vc_expr of the model evaluation)
  clone sizes: N clones, K_eff = 1 / sum p_k^2, CV^2 of clone sizes, largest clone share
  ICC: one-way random-effects ANOVA of the deviations between clones (summed over genes),
       V_clone = (MSB - MSW) / n0, ICC = V_clone / (V_clone + MSW)
  prediction: V_clone x E[sum_k a_k b_k] over the same random splits (a_k, b_k = clone shares of the halves)
Across perturbations: strength against sum p_k^2 (formula: slope = ICC), and measured against predicted clone part.
Usage: choose_clone.py <export dir> <cbc dir> <cbc_libs.tsv> <out.json>"""
import os, sys, glob, gzip, json, collections
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix

MINU, MINC, NSPLIT, SEED = 2, 60, 50, 0


def load_expr(d):
    cells = open(os.path.join(d, "cells.txt")).read().strip().split("\n")
    genes = open(os.path.join(d, "genes.txt")).read().strip().split("\n")
    i = np.fromfile(os.path.join(d, "i.bin"), np.int32); p = np.fromfile(os.path.join(d, "p.bin"), np.int32)
    x = np.fromfile(os.path.join(d, "x.bin"), np.float32)
    X = csr_matrix((x, i, p), shape=(len(cells), len(genes)))          # dgCMatrix columns = cells -> csr rows
    tot = np.asarray(X.sum(1)).ravel()
    X = X.multiply(1e4 / tot[:, None]).tocsr(); X.data = np.log1p(X.data).astype(np.float32)
    meta = pd.read_csv(os.path.join(d, "meta.tsv"), sep="\t", low_memory=False).set_index("cell").loc[cells]
    return cells, np.array(genes), X, meta


def call_clones(cells, cbc_dir, libs_tsv):
    lib_of_run = {}
    t = pd.read_csv(libs_tsv, sep="\t")
    for run, desc in zip(t["Comment[ENA_RUN]"], t["Description"]):
        b, l = desc.replace("CBC, ", "").split(", ")
        if b == "Batch1":
            lib_of_run[run] = int(l.replace("Lib", ""))
    want = {}
    for c in cells:
        cb, lib = c.split("-")[0], int(c.split("_")[-1]); want[(lib, cb)] = c
    reads = collections.Counter()
    for f in glob.glob(os.path.join(cbc_dir, "*.tsv.gz")):
        run = os.path.basename(f).split(".")[0]
        if run not in lib_of_run:
            continue
        lib = lib_of_run[run]
        with gzip.open(f, "rt") as g:
            for line in g:
                _, cb, umi, ucb, n = line.rstrip("\n").split("\t")
                if (lib, cb) in want:
                    reads[(lib, cb, umi, ucb)] += int(n)
    df = pd.DataFrame([(k[0], k[1], k[2], k[3], v) for k, v in reads.items()], columns=["lib", "cb", "umi", "ucb", "reads"])
    stats = {"triples": len(df)}
    df = df[df.reads >= 0.2 * df.groupby(["lib", "cb", "umi"]).reads.transform("max")]
    top = df.groupby(["lib", "cb"]).reads.max().groupby(level=0).median()
    df = df[df.reads >= np.maximum(1, 0.01 * df.lib.map(top))]
    um = df.groupby(["lib", "cb", "ucb"]).size().rename("umis").reset_index()
    # merge UCBs one substitution away from a more abundant one
    tot = um.groupby("ucb").umis.sum().sort_values(ascending=False)
    parent = {}; seen = {}
    for u in tot.index:
        hit = None
        for j in range(len(u)):
            for a in "ACGT":
                if a != u[j]:
                    v = u[:j] + a + u[j + 1:]
                    if v in seen:
                        hit = seen[v]; break
            if hit:
                break
        parent[u] = hit or u
        if not hit:
            seen[u] = u
    um["ucb"] = um.ucb.map(parent); um = um.groupby(["lib", "cb", "ucb"]).umis.sum().reset_index()
    um = um.sort_values(["lib", "cb", "umis"], ascending=[True, True, False])
    first = um.groupby(["lib", "cb"]).head(1).set_index(["lib", "cb"])
    second = um.groupby(["lib", "cb"]).nth(1).set_index(["lib", "cb"]).umis.reindex(first.index).fillna(0)
    ok = (first.umis >= MINU) & (first.umis >= 2 * second)
    clone = {want[k]: first.ucb[k] for k in first.index[ok]}
    stats.update(ucb_before_merge=int(len(tot)), ucb_after_merge=int(len(seen)), cells_with_ucb=int(len(first)),
                 cells_assigned=int(ok.sum()), cells_total=len(cells))
    return clone, stats


def main():
    exp, cbc, libs, outp = sys.argv[1:5]; rng = np.random.default_rng(SEED)
    cells, genes, X, meta = load_expr(exp)
    clone, stats = call_clones(cells, cbc, libs)
    print("clone calls", stats, flush=True)
    meta["clone"] = [clone.get(c) for c in cells]; meta["lib"] = [int(c.split("_")[-1]) for c in cells]
    has = meta.clone.notna().to_numpy()
    # clone-gRNA consistency
    cg = meta[has].groupby("clone").gRNA.agg(["nunique", "size"])
    multi = cg[cg["size"] >= 2]
    stats["clones_2plus_cells"] = int(len(multi)); stats["clones_2plus_one_gRNA"] = float((multi["nunique"] == 1).mean())
    stats["cells_in_multi_gRNA_clones"] = float(multi[multi["nunique"] > 1]["size"].sum() / multi["size"].sum())
    gm = np.asarray(X.mean(0)).ravel() >= 0.1
    Xd = X[:, gm].toarray()
    lib = meta.lib.to_numpy(); P = meta.gRNA.astype(str).to_numpy(); C = meta.clone.to_numpy()
    out = {"rules": __doc__, "clone_calls": stats, "n_genes": int(gm.sum()), "references": {}}
    for refname, grp in (("type", meta.celltype_cl_coarse2.astype(str).to_numpy()), ("library", np.array(["all"] * len(cells)))):
        key = np.array(["%d|%s" % (a, b) for a, b in zip(lib, grp)]); uk, ki = np.unique(key, return_inverse=True)
        S = np.zeros((len(uk), Xd.shape[1])); np.add.at(S, ki, Xd); N = np.bincount(ki, minlength=len(uk)).astype(float)
        V = float(((Xd - S[ki] / N[ki, None]) ** 2).sum() / len(Xd))
        res = {}
        for p in np.unique(P):
            idx = np.where((P == p) & has)[0]
            if len(idx) < MINC:
                continue
            allp = np.where(P == p)[0]                                  # leave the whole perturbation out of its reference
            Sp = np.zeros_like(S); np.add.at(Sp, ki[allp], Xd[allp]); Np = np.bincount(ki[allp], minlength=len(uk))
            ref = (S - Sp) / np.maximum(N - Np, 1)[:, None]
            E = Xd[idx] - ref[ki[idx]]
            cl = C[idx]; uc, ci = np.unique(cl, return_inverse=True); nk = np.bincount(ci).astype(float); pk = nk / nk.sum()
            # random halves and their clone overlap sum_k a_k b_k
            rr, ab = [], []
            for _ in range(NSPLIT):
                bit = rng.integers(0, 2, len(idx)).astype(bool)
                if bit.sum() < 10 or (~bit).sum() < 10:
                    continue
                rr.append(float(E[bit].mean(0) @ E[~bit].mean(0)))
                a = np.bincount(ci[bit], minlength=len(uc)) / bit.sum(); b = np.bincount(ci[~bit], minlength=len(uc)) / (~bit).sum()
                ab.append(float(a @ b))
            cc = []
            for _ in range(NSPLIT):
                side = rng.integers(0, 2, len(uc)).astype(bool)[ci]
                if side.sum() >= 10 and (~side).sum() >= 10:
                    cc.append(float(E[side].mean(0) @ E[~side].mean(0)))
            lu = np.unique(lib[idx]); ll = []
            for _ in range(NSPLIT):
                g = np.isin(lib[idx], rng.permutation(lu)[: len(lu) // 2])
                if g.sum() >= 10 and (~g).sum() >= 10:
                    ll.append(float(E[g].mean(0) @ E[~g].mean(0)))
            # one-way ANOVA between clones
            K = len(uc); Nn = len(idx); m = E.mean(0)
            Mk = np.zeros((K, E.shape[1])); np.add.at(Mk, ci, E); Mk /= nk[:, None]
            ssb = float((nk * ((Mk - m) ** 2).sum(1)).sum()); ssw = float(((E - Mk[ci]) ** 2).sum())
            if K > 1 and Nn > K:
                msb, msw = ssb / (K - 1), ssw / (Nn - K); n0 = (Nn - (nk ** 2).sum() / Nn) / (K - 1)
                vcl = (msb - msw) / n0; icc = vcl / (vcl + msw)
            else:
                vcl = icc = np.nan
            r = dict(n=int(Nn), n_all=int(len(allp)), clones=int(K), K_eff=float(1 / (pk ** 2).sum()), sum_p2=float((pk ** 2).sum()),
                     cv2=float(nk.var() / nk.mean() ** 2), top_clone_share=float(pk.max()), singletons=float((nk == 1).mean()),
                     random=float(np.mean(rr)), clone=float(np.mean(cc)) if cc else np.nan, library=float(np.mean(ll)) if ll else np.nan,
                     ab=float(np.mean(ab)), V_clone=float(vcl), ICC=float(icc), n_libraries=int(len(lu)))
            r["clone_part"] = r["random"] - r["clone"]; r["strength"] = r["clone_part"] / V
            r["predicted_clone_part"] = r["V_clone"] * r["ab"]
            res[p] = r
        tg = {k: v for k, v in res.items() if np.isfinite(v["clone"])}
        a = lambda f: np.array([v[f] for v in tg.values()])
        s, sp2, meas, pred = a("strength"), a("sum_p2"), a("clone_part"), a("predicted_clone_part")
        slope = float((s * sp2).sum() / (sp2 ** 2).sum())
        summ = dict(V_cell=V, n_perturbations=len(tg), mean_random=float(a("random").mean()), mean_clone=float(a("clone").mean()),
                    mean_library=float(np.nanmean(a("library"))), share_clone_part=float(1 - a("clone").mean() / a("random").mean()),
                    random_gt_clone=float((a("random") > a("clone")).mean()), mean_strength=float(s.mean()),
                    strength_se=float(np.std([rng.choice(s, len(s)).mean() for _ in range(2000)])),
                    median_K_eff=float(np.median(a("K_eff"))), median_clones=float(np.median(a("clones"))),
                    median_cv2=float(np.median(a("cv2"))), median_ICC=float(np.nanmedian(a("ICC"))),
                    fit_ICC_slope_strength_on_sum_p2=slope,
                    spearman_strength_sum_p2=float(pd.Series(s).corr(pd.Series(sp2), method="spearman")),
                    pearson_measured_predicted=float(np.corrcoef(meas, pred)[0, 1]),
                    ratio_measured_to_predicted=float(meas.sum() / pred.sum()))
        out["references"][refname] = dict(summary=summ, perturbations=res)
        print(refname, json.dumps(summ), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
