"""CD4 counterpart of the cell-line comparison in cl_ntc_null.py (findings cc): mean random-half cross-half signal of
non-targeting guides divided by that of targeting perturbations, with the same definitions (genes with mean >= 0.1
over non-targeting cells; >= 10 cells per half; targeting effects referenced to all non-targeting cells of the half
and centred over perturbations; non-targeting effects referenced to the mean over kept guides).
Inputs: <ntc_dir> of cd4_ntc_extract.py (sum2.npy, n2.npy) and <cells_dir> of cd4_sc_extract_cells.py (sum8.npy,
n8.npy; index gp * 4 + r1 * 2 + r2; last group = controls). The random bit r1 is the same in both extractions.
Usage: cd4_ntc_ratio.py <out.json> <label>:<ntc_dir>:<cells_dir> [...]"""
import sys, json
import numpy as np


def one(ntc, cells):
    Sn = np.load(ntc + "/sum2.npy"); nn = np.load(ntc + "/n2.npy")
    S8 = np.load(cells + "/sum8.npy", mmap_mode="r"); n8 = np.load(cells + "/n8.npy")
    tot = Sn.sum((0, 1)) / nn.sum()
    genes = json.load(open(ntc + "/meta.json"))["genes"]
    TCR = ("TRAC", "TRBC", "TRAV", "TRBV", "TRGC", "TRDC")
    gk = (tot >= 0.1) & np.array([not g.startswith(TCR) for g in genes])
    ctrl = Sn.sum(0)[:, gk] / nn.sum(0)[:, None]
    okn = (nn >= 10).all(1); Mn = Sn[okn][:, :, gk] / nn[okn][:, :, None]
    Xn = Mn - Mn.mean(0, keepdims=True); rnd = (Xn[:, 0] * Xn[:, 1]).sum(1)
    T = S8.shape[0] - 1                                              # last group = non-targeting controls
    idx = {h: [gp * 4 + h * 2 + r2 for gp in (0, 1) for r2 in (0, 1)] for h in (0, 1)}
    St = np.stack([np.asarray(S8[:T][:, idx[h]][:, :, gk]).sum(1) for h in (0, 1)], 1)
    nt_ = np.stack([n8[:T][:, idx[h]].sum(1) for h in (0, 1)], 1)
    okt = (nt_ >= 10).all(1); Mt = St[okt] / nt_[okt][:, :, None]
    Xt = Mt - ctrl[None]; Xt = Xt - Xt.mean(0, keepdims=True); tc = (Xt[:, 0] * Xt[:, 1]).sum(1)
    return dict(n_genes=int(gk.sum()), n_constructs=int(okn.sum()), random_mean_cross=float(rnd.mean()),
                n_targets=int(okt.sum()), target_random_mean_cross=float(tc.mean()), target_random_median_cross=float(np.median(tc)),
                ratio_nt_to_target_mean=float(rnd.mean() / tc.mean()))


def main():
    out = {"rules": __doc__}
    for arg in sys.argv[2:]:
        lab, ntc, cells = arg.split(":")
        out[lab] = one(ntc, cells); print(lab, json.dumps(out[lab]), flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
