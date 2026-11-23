"""Complete saved zero sets that miss a close pair (PARI lfunzeros at its default subdivision): for every set whose count
drifts from N(t) by more than 1.2, search the largest unfolded gaps (up to 30) again with divz = 256 and add what is found.
usage: fixpairs.py"""
import os, glob, math
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
Z = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')
def drift(g, q):
    g = np.sort(g[g > 20]); u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); r = np.arange(len(g)) - (u - u[0])
    k = max(20, len(g) // 20); return float(r[-k:].mean() - r[:k].mean())
def fix(g, D):
    q = abs(D)
    for _ in range(3):
        if abs(drift(g, q)) <= 1.2: break
        gs = np.sort(g); u = gs / (2 * math.pi) * np.log(q * gs / (2 * math.pi * math.e)); z = []
        for i in np.argsort(-np.diff(u))[:30]:
            z += [float(t) for t in pari('lfunzeros(lfuncreate(%d),[%f,%f],256)' % (D, gs[i] + 1e-6, gs[i + 1] - 1e-6))]
            if z: break
        if not z: break
        g = np.sort(np.concatenate([g, z]))
    return g, drift(g, q)
log = []
for f in ('r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'):
    p = os.path.join(Z, f); z = dict(np.load(p)); changed = False
    for k in list(z):
        D = int(k[1:])
        if abs(drift(z[k], abs(D))) > 1.2:
            n0 = len(z[k]); z[k], d = fix(z[k], D); log.append((f, D, n0, len(z[k]), round(d, 2))); changed = True
    if changed: np.savez_compressed(p, **z)
for p in glob.glob(os.path.join(Z, 'r122', 'D*.npy')) + glob.glob(os.path.join(Z, 'r128', 'D*.npy')):
    D = int(os.path.basename(p)[1:-4]); g = np.load(p)
    if abs(drift(g, abs(D))) > 1.2:
        n0 = len(g); g, d = fix(g, D); np.save(p, g); log.append((p.split('/')[-2], D, n0, len(g), round(d, 2)))
for x in log: print(x)
print(len(log), 'sets fixed or attempted; remaining drift > 1.2:', sum(abs(x[4]) > 1.2 for x in log))
