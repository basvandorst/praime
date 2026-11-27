import numpy as np
from analyze import *
np.seterr(divide='ignore', invalid='ignore')
ns = np.arange(2, 201); kind = np.array([classify(n) for n in ns])
for run in ['runA.npz', 'runB.npz']:
    R = np.load(run); x0 = R['x'][0]; lo, hi = window(x0, 0.2, 0.8)
    for t in [0.2, 0.4, 0.6]:
        it = int(np.argmin(np.abs(R['t'] - t))); g = R['x'][it]
        dn = dirlog(t, 400); l = np.log(ns)
        F = Fsum(g, kaiser(g, lo, hi), ns.astype(float))
        P = template(t, ns, lo, hi, dn)
        seff = np.sum(np.abs(P) ** 2 * np.angle(F / P) * l) / np.sum(np.abs(P) ** 2 * l * l)
        # variant (post hoc): L taken at the starting height of each zero, gamma - s
        gg = np.linspace(lo, hi, 20001); ww = kaiser(gg, lo, hi)
        I2 = np.trapezoid(ww[None, :] * np.exp(-(t / 4) * np.log((gg - seff) / (2 * np.pi))[None, :] * l[:, None]), gg, axis=1)
        P2 = -(1 / (2 * np.pi)) * dn[ns] * ns ** -0.5 * I2
        out = []
        for PP in (P, P2):
            Rr = (F * np.exp(-1j * seff * l) / PP).real
            out.append(' '.join('%s %+.1e' % (k, np.mean(Rr[kind == k]) - 1) for k in ['prime', 'ppow', 'comp']))
        print('%s t=%.1f mean ReR-1 frozen: %s | L at start height: %s' % (run, t, out[0], out[1]))
