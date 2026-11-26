"""Frozen predictions for the blind test (no data used).  For each family (weight k, list of prime levels):
bins of width 0.05 in y = p/N^2 on (0, 3]; for each bin the dim-weighted mean over the (N, p) pairs (p prime <= 3N^2, p != N)
of rho_N (1/sqrt(p) + M(y)), rho_N = (k-1)N/(12 dim), for the formula and for every registered comparator;
predicted bin errors sigma_b = sqrt(c * sum dim p)/sum dim with c = 0.96 (independent-forms model, calibrated on discovery);
power lambda_c = sum_b ((F_formula - F_c)/sigma_b)^2.
Usage: predict.py CODEDIR OUTJSON SETNAME     (SETNAME = holdout or dry)."""
import sys, json, hashlib, datetime
sys.path.insert(0, sys.argv[1])
import numpy as np
from sympy import primerange
from hurwitz import hurwitz12
from trace import Tracer
from model import M_nat, M_parts, M_nat_nu1, M_zub, M_harm
from kernels import K3, K2

SETS = {
    "holdout": [("H2", 2, list(primerange(211, 291))), ("H4", 4, list(primerange(173, 200)))],
    "dry": [("D2", 2, list(primerange(101, 200))), ("D4", 4, list(primerange(101, 170)))],
}
outjson, setname = sys.argv[2], sys.argv[3]
YMAX, BW, CNOISE = 3.0, 0.05, 0.96
H = hurwitz12(100)
grid = np.linspace(1e-5, YMAX + 0.01, 30001)
g3 = np.linspace(1e-4, YMAX + 0.01, 1201)
fams = []
for name, k, levels in SETS[setname]:
    sm, cn = M_parts(grid, k)
    k3 = 0.636446 * np.interp(grid, g3, np.array([K3(z, k) for z in g3]))
    models = {
        "formula": sm + cn,
        "smooth_part": sm,
        "classnumber_part": cn,
        "C1_null": 0 * grid,
        "C2_smooth_only": sm,
        "C3_no_hyperbolic_term": M_nat(grid, k, hyperbolic=False),
        "C4_nu_equal_1": M_nat_nu1(grid, k),
        "C5_harmonic_density": M_harm(grid, k, Dmax=4000),
        "C6_GL3_main_line_only_0.6364K3": k3,
        "C7_degree2_kernel_0.6364K2": 0.636446 * K2(grid, k),
        "C8_Zubrilina_shape": M_zub(grid, k),
    }
    edges = np.arange(0, YMAX + 1e-9, BW); nb = len(edges) - 1
    sw = np.zeros(nb); swp = np.zeros(nb); cnt = np.zeros(nb)
    acc = {m: np.zeros(nb) for m in models}; accb = np.zeros(nb)
    dims = []
    for N in levels:
        dim = Tracer(N, k, H).tr(1); dims.append(dim)
        rho = (k - 1) * N / (12.0 * dim)
        p = np.array([q for q in primerange(2, int(YMAX * N * N) + 1) if q != N], dtype=float)
        y = p / N ** 2; b = np.digitize(y, edges) - 1; ok = (b >= 0) & (b < nb); p, y, b = p[ok], y[ok], b[ok]
        np.add.at(sw, b, dim); np.add.at(swp, b, dim * p); np.add.at(cnt, b, 1)
        np.add.at(accb, b, dim * rho / np.sqrt(p))
        for m, g in models.items():
            np.add.at(acc[m], b, dim * rho * np.interp(y, grid, g))
    base = accb / sw
    pred = {m: (base + acc[m] / sw).tolist() for m in models}
    pred["base_plancherel"] = base.tolist()
    sig = np.sqrt(CNOISE * swp) / sw
    F = np.array(pred["formula"]) - base
    amp_sigma = 1 / np.sqrt(np.sum(F ** 2 / sig ** 2))
    power = {m: float(np.sum(((np.array(pred["formula"]) - np.array(pred[m])) / sig) ** 2)) for m in models if m.startswith("C")}
    S = np.array(pred["smooth_part"]) - base; C = np.array(pred["classnumber_part"]) - base
    X = np.vstack([S, C]).T / sig[:, None]; cov = np.linalg.inv(X.T @ X)
    fams.append(dict(name=name, k=k, levels=levels, dims=dims, bin_edges=edges.tolist(), n_pairs=cnt.astype(int).tolist(),
                     pred=pred, pred_se=sig.tolist(), amp_sigma=float(amp_sigma), power=power,
                     twocomp_sigma=[float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))],
                     twocomp_corr=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]))))
    print(name, "k=%d levels=%d pairs=%d  amp sigma %.3f  two-comp sigmas %.3f %.3f (corr %.2f)" % (
        k, len(levels), int(cnt.sum()), amp_sigma, *fams[-1]["twocomp_sigma"], fams[-1]["twocomp_corr"]))
    for m, v in power.items():
        print("     power vs %-34s lambda = %8.1f" % (m, v))
doc = dict(set=setname, frozen=str(datetime.date.today()), ymax=YMAX, bin_width=BW, noise_factor=CNOISE,
           data_dir="out_holdout" if setname == "holdout" else "out", families=fams,
           decisions=dict(
               P1="k=2: |a-1| <= 3 sigma_a and a/sigma_a >= 8 (a = matched-filter amplitude of the formula, Plancherel term fixed)",
               P2="k=4: same as P1",
               P3="each weight: chi2(formula) over the 60 bins (empirical bin errors) <= 91.95 (99.5% quantile, 60 dof)",
               P4="k=4: two-component fit a*smooth + b*class-number part: |a-1| <= 3 sigma_a, |b-1| <= 3 sigma_b, b/sigma_b >= 5",
               P5="each weight, every comparator with predicted power lambda >= 25: chi2(comparator) - chi2(formula) >= 9",
               S1="each weight: median(empirical/predicted bin error) in [0.85, 1.15] (reported)",
               S2="each weight: amplitudes for N = 1 and N = 3 mod 4 agree within 3 sigma (reported)",
               verdict="CONFIRMED iff P1-P5 all pass; FALSIFIED if a is outside [1-3sigma, 1+3sigma] in either weight, or chi2(formula) > 91.95, or a comparator with lambda >= 25 fits better (delta chi2 < 0)"))
s = json.dumps(doc, indent=1, sort_keys=True)
open(outjson, "w").write(s)
print("sha256", hashlib.sha256(s.encode()).hexdigest(), outjson)
