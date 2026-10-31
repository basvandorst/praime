# Is primality of sieve survivors (sieve to y = x^0.405) a fair independent coin?
# Pair test: P(s_i and s_{i+j} both prime) / rho_emp^2, binned by lag j and by distance d mod 6 classes.
import numpy as np, math
from facts_maxgap_is import PR
def run(X, T=4 * 10**7):
    lo = X; hi = lo + T; y = X ** 0.405
    isp = np.ones(T, bool); surv = np.ones(T, bool)
    for q in PR[PR <= math.isqrt(hi)]:
        q = int(q); st = (-lo) % q; isp[st::q] = False
        if q <= y: surv[st::q] = False
    S = np.nonzero(surv)[0]; b = isp[S].astype(float); rho = b.mean(); n = len(S)
    print(f"x={X:.0e} y={y:.0f}: {n} survivors, rho_emp={rho:.4f}")
    for j in (1, 2, 3, 5, 10):
        both = b[:-j] * b[j:]; r = both.mean() / rho**2
        se = math.sqrt((1 - rho**2) / (rho**2 * (n - j)))
        print(f"  lag {j:2d}: ratio {r:.4f} +- {se:.4f}")
    # by distance d: does the pair factor depend on d (as HL's factors for q > y would)?
    d = S[1:] - S[:-1]; both = b[:-1] * b[1:]
    for dd in (2, 4, 6, 10, 12, 30):
        m = d == dd
        if m.sum() > 2000:
            r = both[m].mean() / (b[:-1][m].mean() * b[1:][m].mean())
            print(f"  next survivor at d={dd:3d}: ratio {r:.4f} (n={m.sum()})")
    # triples: P(three consecutive survivors prime)/rho^3
    t = (b[:-2] * b[1:-1] * b[2:]).mean() / rho**3
    print(f"  triple ratio {t:.4f}")
for X in (10**9, 10**12):
    run(X)
