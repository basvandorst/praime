import numpy as np
from analyze import *
np.seterr(divide='ignore', invalid='ignore')

A = np.load('runA.npz'); Ah = np.load('runA_half.npz'); Bz = np.load('runB.npz')
x0 = A['x'][0]; lo, hi = window(x0, 0.2, 0.8)
ns = np.arange(2, 201); kind = np.array([classify(n) for n in ns])
idxA = (x0 >= lo) & (x0 <= hi)
print('window %.1f..%.1f' % (lo, hi))
for t in [0.2, 0.4, 0.6]:
    it = int(np.argmin(np.abs(A['t'] - t)))
    dn = dirlog(t, 400); P = template(t, ns, lo, hi, dn)
    F = {}
    for name, R in [('A', A), ('Ahalf', Ah), ('B', Bz)]:
        g = R['x'][it]; F[name] = Fsum(g, kaiser(g, lo, hi), ns.astype(float))
    dI = np.abs(F['A'] - F['Ahalf']) / np.abs(P); dT = np.abs(F['A'] - F['B']) / np.abs(P)
    print('t=%.1f integrator |F_A - F_Ahalf|/|P|: max %.1e ; truncation |F_A(M=12000) - F_B(M=30000)|/|P|: max %.1e, rms %.1e'
          % (t, dI.max(), dT.max(), np.sqrt(np.mean(dT ** 2))))
    g = A['x'][it]; gB = Bz['x'][it]
    s = np.mean((g - x0)[idxA]); w0 = kaiser(x0, lo, hi); sw = np.sum(w0 * (g - x0)) / w0.sum()
    sB = np.mean((gB - Bz['x'][0])[(Bz['x'][0] >= lo) & (Bz['x'][0] <= hi)])
    for name in ['A', 'B']:
        # effective shift from the phases of all n <= 200: arg(F/P) = s_eff ln n, weighted by |P|^2
        ph = np.angle(F[name] / P); l = np.log(ns); W = np.abs(P) ** 2
        seff = np.sum(W * ph * l) / np.sum(W * l * l)
        Rr = F[name] * np.exp(-1j * seff * l) / P
        print('   %s: s_eff from phases %.6f (s/t %.6f); plain mean displacement %.6f (A) %.6f (B); taper-weighted %.6f; '
              'after s_eff: ImR rms %.1e, ReR-1 rms prime %.1e comp %.1e'
              % (name, seff, seff / t, s, sB, sw, np.sqrt(np.mean(Rr.imag ** 2)),
                 np.sqrt(np.mean((Rr.real[kind == 'prime'] - 1) ** 2)), np.sqrt(np.mean((Rr.real[kind == 'comp'] - 1) ** 2))))
    # sub-windows (run B, same heights): R per third
    edges = np.linspace(lo, hi, 4)
    print('   sub-windows (run B, thirds; L mid %s): ReR for n=' % ' '.join('%.2f' % np.log(0.5 * (edges[k] + edges[k + 1]) / 2 / np.pi) for k in range(3)))
    sh = [2, 3, 5, 7, 13, 4, 9, 6, 10, 12, 15, 21, 30]
    rows = []
    for k in range(3):
        a, b = edges[k], edges[k + 1]
        Pk = template(t, np.array(sh), a, b, dn)
        Fk = Fsum(gB, kaiser(gB, a, b), np.array(sh, float))
        rows.append(Fk * np.exp(-1j * seff * np.log(sh)) / Pk)
    for j, n in enumerate(sh):
        print('     n=%2d  %s' % (n, '  '.join('%.5f%+.5fi' % (rows[k][j].real, rows[k][j].imag) for k in range(3))))
    # rational (difference) frequencies, relative to the geometric mean of the two parent lines
    rat = [(3, 2), (4, 3), (5, 2), (5, 3), (5, 4), (7, 2), (7, 3), (7, 5)]
    Fr = Fsum(gB, kaiser(gB, lo, hi), np.array([a / b for a, b in rat]))
    print('   rationals |F(a/b)| / sqrt|P(a)P(b)|:', ', '.join('%d/%d %.1e' % (a, b, abs(Fr[i]) / np.sqrt(abs(P[a - 2] * P[b - 2])))
                                                       for i, (a, b) in enumerate(rat)))
