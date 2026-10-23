"""Supplementary Note checks: the propositions behind clone strength, verified symbolically (sympy) and by simulation
recovery (known truth). Model (per perturbation, G genes): cell i of clone k has expression vector
    x_i = mu + c_k + e_i,   c_k ~ N(0, s_c^2 I_G) independent over clones,  e_i ~ N(0, s_e^2 I_G) independent over cells,
and effects are taken against a reference that is fixed (or estimated from many cells). ICC = s_c^2 / (s_c^2 + s_e^2).
P1  For two disjoint sets of cells A, B with clone shares a_k, b_k (sum a_k = sum b_k = 1),
        E[xbar_A . xbar_B] = |mu|^2 + G s_c^2 sum_k a_k b_k
    (cell noise does not enter because A and B share no cell). Independent replicates (no shared clone): |mu|^2.
P2  Random halves: for n_k cells in clone k split by fair coins, E[sum_k a_k b_k] = sum_k p_k^2 (1 + O(1/n)).
    Hence clone strength (random - independent) / (G (s_c^2 + s_e^2)) = ICC sum p_k^2 = ICC / K_eff, K_eff = 1/sum p_k^2,
    independent of the number of cells.
P3  If N0 clones have sizes with mean m and coefficient of variation CV, sum p_k^2 = N0 E[n^2] / (N0 m)^2 (to leading
    order) = (1 + CV^2) / N0, i.e. K_eff = N0 / (1 + CV^2).
P4  Split-half reliability of a null perturbation (correlation over genes between the two half means, then
    Spearman-Brown) tends to 1 as cells grow when s_c > 0 and K is fixed, although its effect in an independent
    experiment is 0: SB extrapolates agreement within the same clones.
P5  Two contexts R, S share clones; clone k shifts R by c_k and S by q c_k. With s = (x_R + x_S)/2 and
    d = (x_R - x_S)/2 from disjoint halves A and B: E[s_A . d_B] = G s_c^2 (1 + q)(1 - q)/4 sum a_k b_k, non-zero for q != 1
    and zero for independent replicates.
Simulation: screens with known s_c, s_e, clone counts and size CV; the estimator used on data (random-half signal minus
the signal between two different null perturbations, divided by the per-cell variance) is compared with ICC sum p_k^2.
Usage: clone_theory_check.py <out.json>"""
import sys, json
import numpy as np
import sympy as sp


def symbolic():
    out = {}
    # P1 with 2 clones, symbolic shares; expectation computed from second moments
    G, sc2, se2, mu2 = sp.symbols("G s_c2 s_e2 mu2", positive=True)
    a1, b1 = sp.symbols("a1 b1", nonnegative=True)
    a = [a1, 1 - a1]; b = [b1, 1 - b1]
    # xbar_A = mu + sum a_k c_k + noise_A ; E[c_k.c_l] = G sc2 delta_kl ; noise terms independent of everything
    E = mu2 + sum(a[k] * b[l] * (G * sc2 if k == l else 0) for k in range(2) for l in range(2))
    target = mu2 + G * sc2 * (a1 * b1 + (1 - a1) * (1 - b1))
    out["P1"] = bool(sp.simplify(E - target) == 0)
    # P2: E[a_k b_k] for n_k cells of clone k among n, split by coins; a_k = X_k / n_A with X_k ~ Bin(n_k, 1/2); to
    # leading order E[a_k b_k] = (n_k/2)(n_k/2) / ((n/2)(n/2)) = p_k^2
    nk, n = sp.symbols("n_k n", positive=True)
    out["P2_leading"] = bool(sp.simplify(((nk / 2) * (nk / 2)) / ((n / 2) * (n / 2)) - (nk / n) ** 2) == 0)
    # P3: sum p^2 = N0 E[n^2] / (N0 m)^2 with E[n^2] = m^2 (1 + CV^2)
    N0, m, cv = sp.symbols("N0 m CV", positive=True)
    out["P3"] = bool(sp.simplify(N0 * m ** 2 * (1 + cv ** 2) / (N0 * m) ** 2 - (1 + cv ** 2) / N0) == 0)
    # P5: s_A contains (1+q)/2 sum a_k c_k, d_B contains (1-q)/2 sum b_k c_k
    q = sp.symbols("q", real=True)
    E5 = sum(((1 + q) / 2 * a[k]) * ((1 - q) / 2 * b[l]) * (G * sc2 if k == l else 0) for k in range(2) for l in range(2))
    out["P5"] = bool(sp.simplify(E5 - G * sc2 * (1 + q) * (1 - q) / 4 * (a1 * b1 + (1 - a1) * (1 - b1))) == 0)
    return out


def sim_screen(rng, n_null, n_cells, n_clones, cv, sc, se, G=300):
    """returns (estimate as on data, truth ICC * mean sum p^2)"""
    shape = 1 / cv ** 2 if cv > 0 else None
    halves, truth = [], []
    for u in range(n_null):
        w = rng.gamma(shape, cv ** 2, n_clones) if cv > 0 else np.ones(n_clones)
        p = w / w.sum(); k = rng.choice(n_clones, n_cells, p=p)
        c = rng.normal(0, sc, (n_clones, G)); X = c[k] + rng.normal(0, se, (n_cells, G))
        pk = np.bincount(k, minlength=n_clones) / n_cells; truth.append((pk ** 2).sum())
        bit = rng.integers(0, 2, n_cells).astype(bool)
        halves.append((X[bit].mean(0), X[~bit].mean(0)))
    rnd = np.mean([a @ b for a, b in halves])
    diff = np.mean([halves[i][0] @ halves[(i + 1) % n_null][1] for i in range(n_null)])
    V = G * (sc ** 2 + se ** 2); icc = sc ** 2 / (sc ** 2 + se ** 2)
    return (rnd - diff) / V, icc * float(np.mean(truth))


def simulate(rng):
    rows = []
    for icc in (0.01, 0.05, 0.2):
        for n_clones, cv in ((5, 0.5), (20, 1.0), (50, 2.0), (200, 1.0)):
            for n_cells in (100, 400):
                sc = np.sqrt(icc); se = np.sqrt(1 - icc)
                est = [sim_screen(rng, 60, n_cells, n_clones, cv, sc, se) for _ in range(3)]
                e = np.mean([x[0] for x in est]); t = np.mean([x[1] for x in est])
                rows.append(dict(ICC=icc, clones=n_clones, CV=cv, cells=n_cells, estimate=float(e), truth=float(t),
                                 K_eff_formula=float(n_clones / (1 + cv ** 2))))
    est = np.array([r["estimate"] for r in rows]); tru = np.array([r["truth"] for r in rows])
    # P4: split-half reliability of a null perturbation with growing cells, fixed clones
    p4 = []
    for n_cells in (50, 200, 800, 3200):
        r = []
        for _ in range(40):
            k = rng.integers(0, 4, n_cells); c = rng.normal(0, np.sqrt(0.05), (4, 300))
            X = c[k] + rng.normal(0, np.sqrt(0.95), (n_cells, 300)); bit = rng.integers(0, 2, n_cells).astype(bool)
            r.append(np.corrcoef(X[bit].mean(0), X[~bit].mean(0))[0, 1])
        rr = float(np.mean(r)); p4.append(dict(cells=n_cells, split_half_r=rr, spearman_brown=2 * rr / (1 + rr)))
    return dict(recovery=rows, recovery_pearson=float(np.corrcoef(est, tru)[0, 1]),
                recovery_slope=float((est * tru).sum() / (tru ** 2).sum()), spearman_brown_vs_cells=p4)


def main():
    rng = np.random.default_rng(7)
    out = {"rules": __doc__, "symbolic": symbolic()}
    out["simulation"] = simulate(rng)
    print(json.dumps(out["symbolic"]), json.dumps({k: v for k, v in out["simulation"].items() if k != "recovery"}), flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
