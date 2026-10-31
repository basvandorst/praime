# Puzzle 1064: heuristic model for the prime factors of G(p) = gcd(2^p + 1, 3^p + 1), calibration and forecast.
# A factor is q = 1 + d p (6 | d) with -2 and -3 both d-th power residues mod q (proof in results/1064.md).
# Model: P(q = 1 + d p divides G(p)) = [local condition on p mod 4] * 4/d^2 * C(d) / ln(1 + d p), where
#   local: 8 | d (all p), or d = 2 (mod 4) and p = d/2 (mod 4); d = 4 (mod 8) never;
#   4/d^2: -2 and -3 each a d-th power with chance 2/d once its quadratic character is forced to +1;
#   C(d): Hardy-Littlewood factor for 1 + d p prime given p prime.
# Pairs d < e on the same p get the joint HL correction; composite G(p) = at least two factors.
# Run: .venv/bin/python results/code/p1064_model.py [--calibrate]   (forecast seconds; calibration about 2 minutes)
import math, sys, numpy as np, gmpy2
from sympy import primerange
SP = list(primerange(3, 20000))
def C(d):
    c = 2.0
    for l in SP: c *= l / (l - 1) if d % l == 0 else l * (l - 2) / (l - 1)**2
    return c
def allowed(d, p4):
    if d % 6: return False
    if d % 8 == 0: return True
    return d % 4 == 2 and p4 == (d // 2) % 4
DM = 1500
Ds = [d for d in range(6, DM + 1, 6) if d % 8 != 4]
Cd = {d: C(d) for d in Ds}
def Pd(d, p): return 4 / d**2 * Cd[d] / math.log(1 + d * p) if allowed(d, p % 4) else 0.0
def joint(d, e):
    j = 1.0
    for l in SP[:400]:
        if d % l and e % l: j *= (l - 1 - (1 if (d - e) % l == 0 else 2)) * (l - 1) / (l - 2)**2
    return j
def pair_weights(p4):
    A = [d for d in Ds if allowed(d, p4)]; W = []
    for i, d in enumerate(A):
        for e in A[i + 1:]:
            w = 16 / (d * d * e * e) * Cd[d] * Cd[e]
            if w > 1e-12: W.append((d, e, w * joint(d, e)))
    return W
W = {1: pair_weights(1), 3: pair_weights(3)}
def expected_composites(X, lo=5):
    u = np.linspace(math.log(lo), math.log(X), 4000); t = np.exp(u); tot = {}
    for p4 in (1, 3):
        f = np.zeros_like(t)
        for d, e, w in W[p4]: f += w / (np.log(1 + d * t) * np.log(1 + e * t))
        tot[p4] = float(np.trapezoid(f * t / np.log(t) / 2, u))
    return tot
if __name__ == "__main__":
    s1, s3 = sum(w for *_, w in W[1]), sum(w for *_, w in W[3])
    print(f"pair weight: p = 1 (mod 4) {s1:.5f}, p = 3 (mod 4) {s3:.5f}; share of p = 3 (mod 4) {s3 / (s1 + s3):.3f}")
    for X in [2243399, 7052113, 10**7, 10**8, 10**9, 10**10]:
        e = expected_composites(X); print(f"X = {X:.3e}: expected composite G(p) for p <= X: {e[1] + e[3]:8.2f}  (p = 3 mod 4: {e[3]:.2f})")
    if "--calibrate" in sys.argv:          # every factor with d <= 300 for p < 10^6, by d class, against the model
        ps = [int(p) for p in primerange(5, 10**6)]; D3 = [d for d in Ds if d <= 300]
        obs = {d: 0 for d in D3}; pred = {d: 0.0 for d in D3}
        for p in ps:
            for d in D3:
                pr = Pd(d, p)
                if pr == 0: continue
                pred[d] += pr; q = 1 + d * p
                if pow(2, p, q) == q - 1 and pow(3, p, q) == q - 1 and gmpy2.is_prime(q): obs[d] += 1
        for d in D3[:8]: print(f"d = {d:3d}: observed {obs[d]:4d}, model {pred[d]:7.1f}")
        print("all d <= 300: observed", sum(obs.values()), "model", round(sum(pred.values()), 1))
