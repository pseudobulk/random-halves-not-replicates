"""Clone ruler (null_ruler.py) on the generality datasets; one loader per dataset.
zhou      Zhou et al. 2023 Nature (GSE216800): OT-I CD8 T cells, transduced in vitro, 4 days, transferred into B16-OVA
          mice, tumour-infiltrating cells at day 7. CellRanger outputs of samples S1-S4 (B1) and S1-S12 (B2); guides by
          feature barcoding. Each vector carries two guides of one gene (dual-guide construct); constructs are found as
          the guide pairs most often detected together. A cell gets the construct with the most guide UMIs if >= 3 UMIs
          and >= 3x the second construct. Unit = construct; null = constructs of non-targeting guides; sample = CellRanger
          sample (reference group and replicate label).
Usage: run_generality.py <dataset> <data dir> <out.json>"""
import os, sys, io, json, tarfile, glob, collections
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy.sparse import vstack, csr_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from null_ruler import lognorm, ruler


def read_10x_tar(path):
    with tarfile.open(path) as t:
        get = lambda n: t.extractfile([m for m in t.getmembers() if m.name.endswith(n)][0]).read()
        feat = pd.read_csv(io.BytesIO(get("features.tsv.gz")), sep="\t", header=None, compression="gzip")
        bc = pd.read_csv(io.BytesIO(get("barcodes.tsv.gz")), header=None, compression="gzip")[0].to_numpy()
        M = mmread(io.BytesIO(__import__("gzip").decompress(get("matrix.mtx.gz")))).T.tocsr()
    return feat, bc, M


def zhou(d):
    tars = sorted(glob.glob(os.path.join(d, "*_S*_B*.tar.gz")))
    parts = []
    for p in tars:
        feat, bc, M = read_10x_tar(p)
        s = os.path.basename(p).split("_", 1)[1].replace(".tar.gz", "")
        parts.append((s, feat, bc, M)); print(s, M.shape, flush=True)
    feat = parts[0][1]; ge = (feat[2] == "Gene Expression").to_numpy(); gd = (feat[2] == "CRISPR Guide Capture").to_numpy()
    gnames = feat[1].to_numpy()[gd]
    for s, f, _, _ in parts:
        assert (f[1].to_numpy() == feat[1].to_numpy()).all(), s
    G = vstack([M[:, gd] for _, _, _, M in parts]).tocsr()
    # constructs: guide pairs most often co-detected (both >= 2 UMIs in a cell)
    co = collections.Counter()
    for i in range(G.shape[0]):
        row = G.getrow(i); j = row.indices[row.data >= 2]
        if len(j) == 2:
            co[tuple(sorted(j))] += 1
    partner = {}
    for (a, b), n in co.most_common():
        if a not in partner and b not in partner:
            partner[a] = b; partner[b] = a
    cons = {}
    for a in range(len(gnames)):
        b = partner.get(a, a); key = "+".join(sorted({gnames[a], gnames[b]})); cons[a] = key
    ck = sorted(set(cons.values())); cidx = {k: i for i, k in enumerate(ck)}
    A = csr_matrix((np.ones(len(gnames)), (np.arange(len(gnames)), [cidx[cons[a]] for a in range(len(gnames))])), shape=(len(gnames), len(ck)))
    C = (G @ A).toarray()
    o = np.sort(C, 1)[:, ::-1]; top = C.argmax(1)
    okc = (o[:, 0] >= 3) & (o[:, 0] >= 3 * o[:, 1])
    unit = np.array([ck[t] if k else None for t, k in zip(top, okc)], dtype=object)
    tgt_of = lambda k: k.split(".")[0]
    target = np.array([tgt_of(u) if u is not None else "" for u in unit])
    is_null = np.array([t == "NTC" for t in target])
    X = vstack([M[:, ge] for _, _, _, M in parts]).tocsr()
    sample = np.concatenate([[s] * len(bc) for s, _, bc, _ in parts])
    info = dict(cells=int(X.shape[0]), assigned=int(okc.sum()), constructs=len(ck),
                paired_guides=int(sum(1 for a in partner)), ntc_constructs=int(sum(1 for k in ck if k.startswith("NTC"))),
                pairs_same_gene=float(np.mean([gnames[a].split(".")[0] == gnames[b].split(".")[0] for a, b in partner.items()])))
    return lognorm(X), sample, unit, target, is_null, sample, info


def knudsen(d, by_lineage=False):
    """Knudsen et al. 2025 Nature (GSE302504): human CAR-T cells of one donor, 14 days ex vivo, then 21 days in NSG
    mice with MM.1S myeloma; bone marrow. X is the authors' log1p-normalised matrix (no raw counts released). Unit =
    single guide (sgRNAs_2+, cells with one guide at >= 2 UMIs); null = intergenic guides; sample = 10x library;
    replicate = mouse (libraries 11_1 and 11_2 are one mouse). by_lineage: reference within sample x CD4 / CD8."""
    import anndata as ad
    a = ad.read_h5ad(os.path.join(d, "GSE302504_adata_106_21_GEX_final.h5ad"))
    g = a.obs["sgRNAs_2+"].astype(str).to_numpy()
    unit = np.array([u if u not in ("none", "multi") else None for u in g], dtype=object)
    target = np.array([u.rsplit("_", 1)[0] if u is not None else "" for u in unit])
    is_null = np.array([t.startswith("INTERGENIC") for t in target])
    mouse = a.obs["mouse"].astype(str).str.split("_").str[0].to_numpy()
    sample = a.obs["sample"].astype(str).to_numpy()
    if by_lineage:
        calls = a.obs["adata_106_21_calls"].astype(str)
        lin = np.where(calls.str.contains("CD4"), "CD4", np.where(calls.str.contains("CD8"), "CD8", "other"))
        sample = np.array([s + "|" + l for s, l in zip(sample, lin)])
    X = a.X.tocsr() if hasattr(a.X, "tocsr") else csr_matrix(a.X)
    info = dict(cells=int(X.shape[0]), assigned=int(sum(u is not None for u in unit)), null_cells=int(is_null.sum()),
                mice=sorted(set(mouse)), reference="sample x lineage" if by_lineage else "sample")
    return X, sample, unit, target, is_null, mouse, info


def fagerberg(d):
    """Fagerberg et al. 2025 (FITdb Zenodo 22212134, Fagerberg_2025.h5ad): P14 CD8 T cells transduced in vitro,
    transferred into mice infected with LCMV Armstrong, spleen at day 8; 6 10x lanes. X = raw counts (dense in the
    file, read in blocks). Unit = true_guide (authors' call; 'inconclusive' and missing dropped); null = NEG-CTRL
    guides; sample = lane (reference group and replicate label; lanes are not known to be separate mice)."""
    import h5py, anndata as ad
    a = ad.read_h5ad(os.path.join(d, "Fagerberg_2025.h5ad"), backed="r")
    obs = a.obs.copy(); a.file.close()
    with h5py.File(os.path.join(d, "Fagerberg_2025.h5ad"), "r") as f:
        D = f["X"]; blocks = []
        for r0 in range(0, D.shape[0], 4000):
            blocks.append(csr_matrix(D[r0:r0 + 4000].astype(np.float32)))
    X = vstack(blocks).tocsr()
    g = obs["true_guide"].astype(str).to_numpy()
    unit = np.array([u if u not in ("nan", "inconclusive") else None for u in g], dtype=object)
    target = np.array([u.rsplit("-", 1)[0] if u is not None else "" for u in unit])
    is_null = np.array([t == "NEG-CTRL" for t in target])
    sample = obs["orig.ident"].astype(str).to_numpy()
    info = dict(cells=int(X.shape[0]), assigned=int(sum(u is not None for u in unit)), null_cells=int(is_null.sum()))
    return lognorm(X), sample, unit, target, is_null, sample, info


def tenx(d, pattern):
    """Generic CellRanger outputs with single-guide feature barcoding (Cao et al. 2026, GSE221578: primary human AML
    cells, transduced, then cultured in vitro or transplanted into mice). A cell gets its top guide if >= 3 UMIs and >= 3x
    the second. Unit = guide; null = guides named neg*, NTC*, ROSA26* / Rosa26* (safe-targeting); sample = file."""
    tars = sorted(glob.glob(os.path.join(d, pattern)))
    Xs, units, samples = [], [], []
    feat0 = None
    for p in tars:
        with tarfile.open(p) as t:
            get = lambda n: t.extractfile([m for m in t.getmembers() if m.name.endswith(n)][0]).read()
            feat = pd.read_csv(io.BytesIO(get("features.tsv.gz")), sep="\t", header=None, compression="gzip")
            M = mmread(io.BytesIO(__import__("gzip").decompress(get("matrix.mtx.gz")))).T.tocsr()
        ge = (feat[2] == "Gene Expression").to_numpy(); gd = (feat[2] == "CRISPR Guide Capture").to_numpy()
        if not gd.any():
            print("no guide features, skipped:", os.path.basename(p), flush=True); continue
        if feat0 is None:
            feat0 = feat[1].to_numpy()[ge]
        assert (feat[1].to_numpy()[ge] == feat0).all()
        G = M[:, gd].toarray(); names = feat[1].to_numpy()[gd]
        o = np.sort(G, 1)[:, ::-1]; okc = (o[:, 0] >= 3) & (o[:, 0] >= 3 * o[:, 1])
        units.append(np.array([names[i] if k else None for i, k in zip(G.argmax(1), okc)], dtype=object))
        Xs.append(M[:, ge]); s = os.path.basename(p).split("_filtered")[0]; samples.append(np.array([s] * M.shape[0]))
        print(s, M.shape, int(okc.sum()), flush=True)
    unit = np.concatenate(units); sample = np.concatenate(samples)
    isn = lambda u: u is not None and (u.lower().startswith("neg") or u.startswith("NTC") or u.lower().startswith("rosa26"))
    is_null = np.array([isn(u) for u in unit])
    target = np.array(["NULL" if isn(u) else (u.split("_")[0].split(".")[0] if u is not None else "") for u in unit])
    X = vstack(Xs).tocsr()
    info = dict(cells=int(X.shape[0]), assigned=int(sum(u is not None for u in unit)), null_cells=int(is_null.sum()),
                samples=sorted(set(sample)))
    return lognorm(X), sample, unit, target, is_null, None, info


def schmidt(d, cond):
    """Schmidt et al. 2022 Science (FITdb Zenodo 22212134, Schmidt_2022_resting / _restim.h5ad): primary human T cells
    of two donors mixed 1:1, CRISPRa library transduced together (MOI 0.3), then rested or restimulated; donors assigned
    by genotype. X = raw counts. Unit = guide_id; null = NO-TARGET guides; sample = donor (reference group and replicate:
    the two donors share no clones)."""
    import anndata as ad
    a = ad.read_h5ad(os.path.join(d, "Schmidt_2022_%s.h5ad" % cond))
    keep = (a.obs["donor"].astype(str) != "Unassigned").to_numpy(); a = a[keep]
    X = a.X.tocsr() if hasattr(a.X, "tocsr") else csr_matrix(a.X)
    unit = a.obs["guide_id"].astype(str).to_numpy().astype(object)
    target = a.obs["gene"].astype(str).to_numpy(); is_null = target == "NO-TARGET"
    donor = a.obs["donor"].astype(str).to_numpy()
    info = dict(cells=int(X.shape[0]), null_cells=int(is_null.sum()), null_guides=int(len(set(unit[is_null]))), condition=cond)
    return lognorm(X.astype(np.float32)), donor, unit, target, is_null, donor, info


def main():
    ds, d, outp = sys.argv[1:4]
    loaders = {"zhou": zhou, "knudsen": knudsen, "knudsen_lineage": lambda d: knudsen(d, True), "fagerberg": fagerberg,
               "schmidt_resting": lambda d: schmidt(d, "resting"), "schmidt_restim": lambda d: schmidt(d, "restim"),
               "cao_invivo": lambda d: tenx(d, "*invivo*.tar.gz"), "cao_invitro": lambda d: tenx(d, "*invitro*.tar.gz")}
    X, sample, unit, target, is_null, rep, info = loaders[ds](d)
    print(ds, info, flush=True)
    out = dict(rules=__doc__, dataset=ds, info=info)
    for mode in ("all", "null"):
        out["ruler_" + mode] = ruler(X, sample, unit, target, is_null, rep, ref_mode=mode)
        print(mode, json.dumps(out["ruler_" + mode]), flush=True)
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
