"""Loop 6: high zeros of zeta (Odlyzko, near zero 10^12, 10^21, 10^22) moved by the universal flow, against flowed Haar matrices."""
import sys, numpy as np
from scipy.stats import unitary_group
from loop3 import flowed_angles
BASE = {3: 267653395647.0, 4: 144176897509546973000.0, 5: 1370919909931995300000.0}
def load(f):
    lines = open('/home/bas/praime/writeup/data/riemann/odlyzko_zeros%d.txt' % f).read().split('\n')
    v = []
    for ln in lines:
        s = ln.strip()
        try: v.append(float(s))
        except ValueError: pass
    v = np.array(v); v = v[np.concatenate(([True], np.diff(v) > 0))] if False else v
    return v[-10000:]
def rhs(x, a, b):
    out = np.empty_like(x)
    for i0 in range(0, len(x), 1000):
        d = x[i0:i0 + 1000, None] - x[None, :]
        np.fill_diagonal(d[:, i0:i0 + 1000], np.inf)
        out[i0:i0 + 1000] = np.sum(1.0 / d, axis=1)
    return (out + np.log((b - x) / (x - a))) / (2 * np.pi ** 2)
def flow(x, taus):
    a, b = x[0] - 0.5, x[-1] + 0.5; tau = 0.0; snaps = []; n = 0
    for T in taus:
        while tau < T - 1e-12:
            g = np.min(np.diff(x)); dt = min(0.08 * np.pi ** 2 * g * g, 0.02, T - tau)
            k1 = rhs(x, a, b); k2 = rhs(x + dt * k1, a, b)
            x = x + 0.5 * dt * (k1 + k2); tau += dt; n += 1
        snaps.append(x.copy())
    return snaps, n
if __name__ == "__main__":
    taus = [0.0, 0.5, 1.0, 2.0, 4.0]
    res = {T: [] for T in taus}; amp = {T: [] for T in taus}
    for f in (3, 4, 5):
        g = load(f); L = np.log(BASE[f] / (2 * np.pi)); x_all = g * L / (2 * np.pi)
        for blk in range(3):
            x0 = x_all[blk * 3333:(blk + 1) * 3333]
            snaps, n = flow(x0.copy(), taus)
            for T, x in zip(taus, snaps):
                mid = x[667:2667]; s = np.diff(mid)
                res[T].append(s.var() / s.mean() ** 2)
                amp[T].append([abs(np.sum(np.exp(2j * np.pi * mid * np.log(p) / L))) for p in (2, 3, 5)])
            print(f"file {f} block {blk}: {n} steps", flush=True)
    N = 200; zs = [np.linalg.eigvals(U) for U in unitary_group.rvs(N, size=150, random_state=np.random.default_rng(21))]
    print("  tau   zeta: spacing variance (9 blocks)   matrices N = 200     z      line amplitude ratio to tau = 0, p = 2, 3, 5 (predicted)")
    a0 = np.array(amp[0.0])
    for T in taus:
        v = np.array(res[T]); m, e = v.mean(), v.std(ddof=1) / 3
        sp = []
        for z in zs:
            ang, _ = flowed_angles(z, T / N ** 2)
            if len(ang) == N: sp.append(np.diff(np.append(ang, ang[0] + 2 * np.pi)) * N / (2 * np.pi))
        sp = np.array(sp); mv = sp.var(); me = sp.var(axis=1).std(ddof=1) / np.sqrt(len(sp))
        r = np.array(amp[T]) / a0
        Ls = np.repeat([np.log(BASE[f] / (2 * np.pi)) for f in (3, 4, 5)], 3)
        pred = [np.mean(np.exp(-T * (np.log(p) / Ls) * (1 - np.log(p) / Ls))) for p in (2, 3, 5)]
        print(f" {T:4.1f}   {m:7.4f} +- {e:6.4f}                  {mv:7.4f} +- {me:6.4f}  {(m-mv)/np.hypot(e,me):+5.1f}   "
              + "  ".join(f"{r[:, i].mean():6.4f}+-{r[:, i].std(ddof=1)/3:6.4f} ({pred[i]:6.4f})" for i in range(3)))
