"""Loop 8: binned form factor of the flowed high zeros of zeta against the cluster series K_tau(u)."""
import numpy as np
from loop6 import load, flow, BASE
from loop2 import K_limit
taus = [0.0, 1.0, 4.0]; edges = np.arange(0.05, 1.56, 0.1); ug = np.arange(0.05, 1.55, 0.0005)
out = {T: [] for T in taus}
for f in (3, 4, 5):
    g = load(f); L = np.log(BASE[f] / (2 * np.pi)); x_all = g * L / (2 * np.pi)
    for blk in range(3):
        snaps, n = flow(x_all[blk * 3333:(blk + 1) * 3333].copy(), taus)
        for T, x in zip(taus, snaps):
            mid = x[667:2667]; mid = (mid - mid[0]) / np.mean(np.diff(mid))
            K = np.array([np.abs(np.exp(2j * np.pi * ug[i:i + 200, None] * mid[None, :]).sum(axis=1)) ** 2 for i in range(0, len(ug), 200)]).ravel() / len(mid)
            out[T].append([K[(ug >= a) & (ug < b)].mean() for a, b in zip(edges[:-1], edges[1:])])
print("  u bin        " + "".join(f"tau = {T:3.1f}: zeta (series)        " for T in taus))
for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
    row = f" {a:4.2f}-{b:4.2f}  "
    for T in taus:
        v = np.array(out[T])[:, i]; m, e = v.mean(), v.std(ddof=1) / 3
        if b <= 1.0001:
            us = np.linspace(a, b, 5)[:-1] + (b - a) / 8
            s = np.mean([sum(K_limit(u, T, rmax=5, n=60000)) for u in us])
            row += f"  {m:7.4f}+-{e:6.4f} ({s:7.4f}) z={(m-s)/e:+4.1f}"
        else:
            row += f"  {m:7.4f}+-{e:6.4f}                 "
    print(row)
