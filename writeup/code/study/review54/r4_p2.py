# Prime-square secondary term for level-1 forms (Q = 1): e(-2 sqrt(p^2)) = 1 makes n = p^2 coherent.
# I_n = -(1/2pi) Lambda(n) n^{-1/2} int_0^T exp(i(2 theta(t) - t ln n)) dt   (explicit formula, sharp window)
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, numpy as np, cypari2
from r1_sums import theta, own_data, D
pari = cypari2.Pari(); pari.allocatemem(200*10**6)
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
L = own_data()
import sympy
out = {}
for key in ('D12','F16','F18','F20','F22','F26'):
    v = L[key]; P = pred[key]; T2 = float(P['T2']); k = int(round(v['mu'][0]*2+1))
    co = [int(x) for x in pari('mfcoefs(mfbasis(mfinit([1,%d],1))[1], 400)' % k)]
    assert co[1] == 1
    lam = lambda p: co[p] / p ** ((k - 1) / 2)
    t = np.linspace(1e-6, T2, 600001); dt = t[1] - t[0]
    th2 = 2 * theta(t, v['Q'], v['mu'])
    tot = {'p2': 0j, 'p2_mean_only': 0j, 'p3': 0j}
    contrib = []
    for p in sympy.primerange(2, 400):
        lp = lam(p)
        for m, name in ((2, 'p2'), (3, 'p3')):
            n = p ** m
            if np.log(n) > 2 * 8 + 20: continue
            # Lambda(p^m) = (alpha^m + beta^m) log p ; alpha+beta = lam, alpha beta = 1
            pm = {2: lp * lp - 2, 3: lp ** 3 - 3 * lp}[m]
            f = np.exp(1j * (th2 - t * np.log(n)))
            I = (f.sum() - 0.5 * (f[0] + f[-1])) * dt
            val = -(1 / (2 * np.pi)) * pm * np.log(p) * n ** -0.5 * I
            tot[name] += val
            if m == 2:
                tot['p2_mean_only'] += -(1 / (2 * np.pi)) * (-1) * np.log(p) * n ** -0.5 * I
    eps = v['eps']; c = P['c']; pr = -eps * c * T2 / (2 * np.pi)
    g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']; g = np.sort(g[(g > 0) & (g <= T2)])
    G = np.exp(2j * theta(g, v['Q'], v['mu'])).sum()
    rot = np.conj(-eps)
    out[key] = dict(G=[G.real, G.imag], pred=[pr.real, pr.imag], P2=[tot['p2'].real, tot['p2'].imag], P3=[tot['p3'].real, tot['p3'].imag], P2mean=[tot['p2_mean_only'].real, tot['p2_mean_only'].imag])
    print('%s k=%d  G-pred (rotated) %+.1f%+.1fi | P2 rotated %+.1f%+.1fi (mean-only %+.1f%+.1fi) | P3 rotated %+.1f%+.1fi | slope %.3f -> corrected %.3f'
          % (key, k, ((G - pr) * rot).real, ((G - pr) * rot).imag, (tot['p2'] * rot).real, (tot['p2'] * rot).imag,
             (tot['p2_mean_only'] * rot).real, (tot['p2_mean_only'] * rot).imag, (tot['p3'] * rot).real, (tot['p3'] * rot).imag,
             (G / pr).real, ((G - tot['p2']) / pr).real))
json.dump(out, open(os.path.dirname(os.path.abspath(__file__))+'/r4_p2.json', 'w'), indent=1)
