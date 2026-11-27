"""Loop 5: spacing variance of the flowed zeros of zeta against flowed Haar matrices at the same tau = t L^2 / 4."""
import sys, numpy as np
from scipy.stats import unitary_group
from loop3 import flowed_angles
d = np.load(sys.argv[1]); ts = d['t']; X = d['x']; lo, hi = 6000, 10000
N = 100; zs = [np.linalg.eigvals(U) for U in unitary_group.rvs(N, size=400, random_state=np.random.default_rng(11))]
def mat_var(tau):
    sp = []
    for z in zs:
        a, _ = flowed_angles(z, tau / N ** 2)
        if len(a) == N: sp.append(np.diff(np.append(a, a[0] + 2 * np.pi)) * N / (2 * np.pi))
    sp = np.array(sp); return sp.var(), sp.var() * np.sqrt(2.0 / sp.size) * 1.5, np.mean((sp - 1) ** 3)
print("   t     tau    zeta: var   skew3     matrices: var      skew3    ratio")
for ti in (0, 1, 2, 4, 6, 9, 12):
    g = X[ti, lo:hi]; L = np.log(g / (2 * np.pi)); s = np.diff(g) * L[:-1] / (2 * np.pi); s = s / s.mean()
    tau = ts[ti] * np.mean(L ** 2) / 4
    v, e, k3 = mat_var(tau)
    print(f" {ts[ti]:5.2f} {tau:6.2f}   {s.var():7.4f}  {np.mean((s-1)**3):+7.4f}    {v:7.4f}+-{e:6.4f}  {k3:+7.4f}   {s.var()/v:5.3f}")
