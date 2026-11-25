"""freeze.py -- write predictions_twinzeros.json (predictions, alternatives, decision rules, expected errors, hashes)."""
import json, hashlib, os
import theory

HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = "/home/bas/praime/writeup/data/riemann/"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

dry = json.load(open(os.path.join(HERE, "tz_test_result_dry.json")))["stats"]
S2 = theory.singular_const("h2")
# expected holdout errors: dry-run errors on 2e8..2e9 (same length in ln x) times the ratio of the mean ln x
# (noise of a k-term pattern in Y grows like (ln x)^(k/2)): pairs 22.57/20.27 = 1.11, triples 1.11^1.5 = 1.17
f2, f3 = 22.57 / 20.27, (22.57 / 20.27) ** 1.5
exp_err = {}
for key, v in dry.items():
    fac = f3 if key.startswith("T0") else f2
    if isinstance(v, dict):
        for lab, x in v.items():
            if isinstance(x, list) and len(x) == 3:
                exp_err["%s | %s" % (key, lab)] = {"predicted": x[2], "expected_sigma": round(x[1] * fac, 4)}
    else:
        exp_err["%s | raw chi_-3" % key] = {"predicted": v[2], "expected_sigma": round(v[1] * fac, 4)}

pred = {
    "title": "Linearised explicit formula for prime patterns: frozen holdout test",
    "frozen_on": "2026-11-25",
    "discovery_range": "x <= 2e9 (pairs, triples, Sophie Germain pairs); prime sums (regressors) to 4e9",
    "holdout_range": "x in (2e9, 2e10]; pair statistics computed on (2e9, 2e10], prime sums on (4e9, 4e10]",
    "statistic": ("Y(u) = (sum_{p<=x, pattern} prod ln(L_i(p)) - S x)/sqrt(x) on the grid u = ln X_k, X_k = floor(exp(k/2048)); "
                  "regressors (theta(a x) - a x)/sqrt(x) and theta(a x, chi)/sqrt(x) for primitive chi of conductor <= 30; "
                  "fit: differences of the series weighted by (ln x)^(-k/2), FFT (rectangular window), bins with 2 <= gamma <= 6000 "
                  "(and 2..300), weighted least squares on real and imaginary parts with nuisance regressors (u, u^2, exp(-u/2), exp(-u), "
                  "constant drift); noise variance per bin from the smoothed squared residuals (41 bins, 3 passes); errors from (X^T W X)^-1. "
                  "Model 'all': zeta coefficient(s) free, one loudness factor for all predicted characters together, free coefficients for the "
                  "null set 3.1, 4.1, 5.2, 7.1, 7.4, 8.1, 8.2, 11.1, 12.1, 9.1 (those not predicted present). Model 'cond': one loudness per conductor."),
    "conjecture": ("sum_{n<=x} prod_i Lambda(a_i n + b_i) = S x + S sum_i (1/a_i)[(psi(a_i x) - a_i x) + sum_{chi primitive} c_i(chi) psi(a_i x, chi)] "
                   "+ (terms not linear in the zeros), c_i(chi) = mean of conj(chi(a_i n + b_i)) over the admissible n mod cond(chi). "
                   "Pairs (n, n+h): zeta coefficient 2 S(h); chi enters with 2 S(h) mu(q) conj(chi(h)) / prod_{p|q}(p-2) iff chi is even with "
                   "odd squarefree conductor q coprime to h, else 0."),
    "predictions_and_expected_errors": exp_err,
    "alternatives": {
        "zeta with factor 1 (coefficient S(h) instead of 2 S(h))": {"h2 zeta": S2},
        "zeta absent": {"h2 zeta": 0.0},
        "no L-function lines": {"loudness": 0.0},
        "proposer's original guess: chi_-3 (mod 3) heard by twins at the loudness of zeta": {"h2 null 3.1": "+-2S = +-%.4f" % (2 * S2)},
    },
    "decisions": {
        "P1": "twins h=2, band [2,6000], model all: |beta_zeta - 2S(2)| <= 3 sigma AND |beta_zeta - S(2)| > 5 sigma",
        "P2": "twins: |loudness - 1| <= 3 sigma AND loudness/sigma >= 5",
        "P3": "twins: |beta(chi_-3)| <= 3 sigma AND min |beta(chi_-3) -+ 2S(2)| > 5 sigma",
        "P4": "twins: chi^2 of the 13 null coefficients (with their covariance) <= 34.53 (99.9% quantile, 13 dof)",
        "S1": "gaps 4, 6, 8, 10, 12, 30: zeta coefficient and loudness, all 12 |z| <= 3.5",
        "S2": "Sophie Germain (n, 2n+1): zeta(a=1) = S, zeta(a=2) = S/2, loudness = 1, chi q=3 and q=4 loudness = 1: all |z| <= 3 and chi q=3 loudness/sigma >= 5",
        "S3": "triples (0,2,6), (0,4,6): chi q=3 loudness |z| <= 3 and the raw chi_-3 coefficient has the predicted sign (-S3, +S3) at >= 2 sigma, both",
        "S4": "twins, band [2,300]: zeta and loudness within 3 sigma (low power)",
        "S5": "twins on the explicit zeros (combs of zeta, chi_5, chi_-3, chi_-4 zeros up to 2500), band [2,1900]: predictions 2S, 2S/3, 0, 0, all |z| <= 3",
        "S6": "twins per conductor 5, 7, 15, 21: loudness 1 within 3 sigma",
        "verdict": "PRIMARY CONFIRMED iff P1-P4 all pass; S1-S6 are reported separately",
    },
    "falsifier": ("If on (2e9, 2e10] the twin (h=2) zeta coefficient misses 2S(2) = 2.6406 by more than 3 sigma, or the joint "
                  "L-function loudness misses 1 by more than 3 sigma, or the mod-3 character coefficient of twins exceeds 3 sigma, "
                  "the linearised explicit formula for prime pairs is falsified."),
    "power": ("expected sigma (holdout): zeta(h2) ~%.3f, loudness ~%.3f, chi_-3 ~%.3f; separation from 'factor 1' ~%.1f sigma, "
              "from 'no L-function lines' ~%.1f sigma, from the chi_-3 = +-2S guess ~%.1f sigma; S3 sign tests have power about 0.85 each" % (
                  exp_err["h2/full/all | zeta(a=1)"]["expected_sigma"], exp_err["h2/full/all | all present chi (loudness)"]["expected_sigma"],
                  exp_err["h2/full/all | null 3.1"]["expected_sigma"],
                  S2 / exp_err["h2/full/all | zeta(a=1)"]["expected_sigma"],
                  1 / exp_err["h2/full/all | all present chi (loudness)"]["expected_sigma"],
                  2 * S2 / exp_err["h2/full/all | null 3.1"]["expected_sigma"])),
    "hashes": {n: sha(os.path.join(HERE, n)) for n in ("tz_counts.py", "sieve_core.py", "tzlib.py", "models.py", "theory.py", "disc.npz")},
}
for z in ("zeros_40000_flint.txt", "lzeros/D5.txt", "lzeros/D-3.txt", "lzeros/D-4.txt"):
    pred["hashes"][ZDIR + z] = sha(ZDIR + z)
json.dump(pred, open(os.path.join(HERE, "predictions_twinzeros.json"), "w"), indent=1)
h = sha(os.path.join(HERE, "predictions_twinzeros.json"))
open(os.path.join(HERE, "predictions_twinzeros.json.sha256"), "w").write(h + "  predictions_twinzeros.json\n")
print(h)
print(pred["power"])
