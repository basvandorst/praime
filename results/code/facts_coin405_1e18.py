import numpy as np, math, gmpy2
from facts_maxgap_is import PR, logprod, omega, GAMMA
X = 10**18; T = 6 * 10**6; y = X ** 0.405
rho_th = 1.0 / (math.log(X) * math.exp(logprod(1, y)) * math.exp(GAMMA) * omega(1 / 0.405))
allb = []
for w in range(4):
    lo = X + w * 10**15
    surv = np.ones(T, bool)
    for q in PR[PR <= y]:
        q = int(q); surv[(-lo) % q::q] = False
    S = np.nonzero(surv)[0]
    allb.append(np.array([gmpy2.is_bpsw_prp(lo + int(s)) for s in S], float))
b = np.concatenate(allb); rho = b.mean(); n = len(b)
print(f"x=1e18, y={y:.3g}: {n} survivors, rho_emp={rho:.4f}, rho_theory={rho_th:.4f}")
for j in (1, 2, 3, 5):
    pr = np.concatenate([bb[:-j] * bb[j:] for bb in allb]).mean() / rho**2
    print(f"  lag {j}: ratio {pr:.4f} +- {math.sqrt((1-rho**2)/(rho**2*n)):.4f}")
t = np.concatenate([bb[:-2]*bb[1:-1]*bb[2:] for bb in allb]).mean() / rho**3
print(f"  triple ratio {t:.4f}")
