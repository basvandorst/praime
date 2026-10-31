# LO-S last-digit transitions: real vs the c = 0.405 rule (next prime = J-th survivor of sieve to y = x^0.405).
import numpy as np, math
from facts_maxgap_is import PR, logprod, omega, GAMMA
def run(X, T=10**7):
    lo = X; L = T; hi = lo + L
    isp = np.ones(L, bool); surv = np.ones(L, bool); y = X ** 0.405
    for q in PR[PR <= math.isqrt(hi)]:
        q = int(q); st = (-lo) % q
        isp[st::q] = False
        if q <= y: surv[st::q] = False
    rho = 1.0 / (math.log(X) * math.exp(logprod(1, y)) * math.exp(GAMMA) * omega(1 / 0.405))
    P = lo + np.nonzero(isp)[0]; S = lo + np.nonzero(surv)[0]
    real = np.zeros((10, 10)); np.add.at(real, (P[:-1] % 10, P[1:] % 10), 1)
    idx = np.searchsorted(S, P[:-200])            # every prime is a survivor (all q <= y < p)
    assert np.all(S[idx] == P[:-200])
    mod = np.zeros((10, 10))
    for j in range(1, 80):
        w = rho * (1 - rho) ** (j - 1)
        np.add.at(mod, (P[:-200] % 10, S[idx + j] % 10), w)
    real /= real.sum(); mod /= mod.sum()
    d = [1, 3, 7, 9]
    print(f"x = {X:.0e}, y = {y:.0f}, rho = {rho:.4f}, {len(P)} primes")
    print("  from->to   real     model    ratio")
    worst = 0
    for a in d:
        for b in d:
            r = mod[a, b] / real[a, b]; worst = max(worst, abs(r - 1))
            print(f"   {a}->{b}    {real[a,b]:.5f}  {mod[a,b]:.5f}  {r:.4f}")
    diag_r = sum(real[a, a] for a in d); diag_m = sum(mod[a, a] for a in d)
    print(f"  diagonal share real {diag_r:.4f} model {diag_m:.4f}; worst cell {100*worst:.2f} %")
for X in (10**9, 10**12):
    run(X)
