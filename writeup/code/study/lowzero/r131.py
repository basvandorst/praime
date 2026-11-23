"""Round 131 (preregister/round131.txt): races b vs 1 mod prime q from the zeros of all characters mod q. usage: r131.py"""
import os, json, math, cmath
import numpy as np, cypari2
from sympy import primerange
from scipy.stats import pearsonr
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
P = np.array(list(primerange(2, 10 ** 6)), dtype=np.int64); lnP = np.log(P.astype(float))
out = []
for q in (5, 7, 11, 13, 17, 19, 23, 29, 31):
    g0 = int(pari('znprimroot(%d)' % q).lift()); ind = {}; v = 1
    for e in range(q - 1): ind[v] = e; v = v * g0 % q
    chars = []
    for k in range(1, q - 1):                                   # chi_k(g0^e) = exp(2 pi i k e/(q-1))
        L = pari.lfuncreate(pari('[znstar(%d,1),[%d]]' % (q, k)))
        z = np.array([float(t) for t in pari.lfunzeros(L, pari('[-50,50]'))])
        real = (2 * k) % (q - 1) == 0
        S = np.zeros(len(P), complex)
        for i0 in range(0, len(z), 40):
            zz = z[i0:i0 + 40]; S += (np.exp(1j * np.outer(lnP, zz)) / (0.5 + 1j * zz)[None, :]).sum(axis=1)
        chars.append((k, real, -S - (1.0 if real else 0.0)))
    # measured race: pi(p;q,b) - pi(p;q,1) after each prime p
    cls = P % q
    for b in range(2, q):
        chi_b = lambda k: cmath.exp(-2j * math.pi * k * ind[b] / (q - 1))
        E = np.zeros(len(P))
        for k, real, T in chars:
            E += np.real((chi_b(k) - 1) * T)
        E /= (q - 1)
        lead = np.cumsum((cls == b).astype(int) - (cls == 1).astype(int))
        m = P > q
        out.append({'q': q, 'b': b, 'pred': float(np.mean(E[m] > 0)), 'meas': float(np.mean(lead[m] > 0))})
    print(q, [(r['b'], round(r['pred'], 2), round(r['meas'], 2)) for r in out if r['q'] == q], flush=True)
pr = np.array([r['pred'] for r in out]); me = np.array([r['meas'] for r in out])
r = pearsonr(pr, me)[0]; mae = float(np.mean(np.abs(pr - me)))
print('races %d; Pearson %.3f; mean abs error %.3f; readings %s %s' % (len(out), r, mae, r >= 0.9, mae <= 0.08))
json.dump({'races': out, 'pearson': float(r), 'mae': mae, 'readings': [bool(r >= 0.9), bool(mae <= 0.08)]}, open(os.path.join(H, 'r131.json'), 'w'), indent=0)
