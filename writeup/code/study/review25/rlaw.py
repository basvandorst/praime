"""Runs r89's per-d law (code copied from r89.law, stored per lam) and compares the z^0 coefficient with r59.share
and with r59.json's stored law. Saves out[k][lam][d] to rlaw.npz."""
import sys, os, json, math, time
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np
import r89, r59
t0 = time.time()
LAMS = r89.LAMS; Z = r89.Z; NZ = r89.NZ
dl = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else list(range(6, 1201, 6))
out = np.zeros((4, 3, len(dl)))
lnS02 = sum(math.log(1 - 2 / q) - 2 * math.log(1 - 1 / q) for q in r59.MID) + r59.TAIL[2]
r59j = json.load(open('/home/bas/praime/writeup/data/study/r59.json'))['law']
worst = 0; worst2 = 0
for di, d in enumerate(dl):
    h, lnB, lns, lnr = r59.lnS_parts(d)
    cls = r59.CLS[(np.gcd(r59.CLS + d, r59.PP) == 1) & (np.gcd(r59.CLS + d + 2, r59.PP) == 1)]
    if len(cls) == 0:
        continue
    M = ((np.gcd(cls[:, None] + h[None, :], r59.PP) == 1) & (np.gcd(cls[:, None] + h[None, :] + 2, r59.PP) == 1)).astype(float) if len(h) else np.zeros((len(cls), 0))
    R = np.expm1(lnr)
    for li, lam in enumerate(LAMS):
        rho0 = 1 / (lam * r59.DEN); p = rho0 ** 2 * np.exp(lns)
        logG = np.zeros((len(cls), NZ), complex)
        if len(h):
            logG += M @ np.log1p(p[:, None] * (Z[None, :] - 1))
        if len(h) > 1:
            W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))
            for m in range(NZ):
                X = M * W[:, m][None, :]
                logG[:, m] += 0.5 * np.sum((X.real @ R + 1j * (X.imag @ R)) * X, axis=1)
        coef = np.fft.fft(np.exp(logG), axis=1).real / NZ
        for k in range(4):
            out[k, li, di] = rho0 ** 2 * math.exp(lnB - lnS02) * coef[:, k].sum() / len(r59.CLS)
    ref = r59j[str(d)]
    rel = max(abs(out[0, li, di] / ref[li] - 1) for li in range(3))
    worst = max(worst, rel)
    if d in (6, 60, 300, 600, 1200) or len(dl) < 20:
        sh = r59.share(d, LAMS)
        rel2 = max(abs(out[0, li, di] / sh[li] - 1) for li in range(3)); worst2 = max(worst2, rel2)
        print(d, 'z0/share-1 %.2e' % rel2, 'z0/r59json-1 %.2e' % rel, [round(x, 10) for x in out[:, 1, di]], round(time.time() - t0), 's', flush=True)
print('max rel diff z0 vs r59.json', worst, 'vs share', worst2)
np.savez('rlaw.npz' if len(dl) > 20 else 'rlaw_small.npz', out=out, d=np.array(dl), lams=np.array(LAMS))
