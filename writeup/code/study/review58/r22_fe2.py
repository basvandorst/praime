"""FE/theta check e^{2i theta(t)} L(1/2+it)/conj L(1/2+it) = eps for curves with 2 | N or a non-squarefree conductor."""
import sys, numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import theta, load
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
rows = {r['key']: r for r in load() if r['cls'] == 'g2'}
for k in ('G2_N1532_38', 'G2_N2862_17', 'G2_N2758_36', 'G2_N2575_6'):
    r = rows[k]; L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2])); Li = pari.lfuninit(L, [1, 100])
    out = []
    for t in (3.3, 41.7, 97.9):
        v = complex(pari.lfun(Li, pari('1 + %s*I' % t)))
        out.append(np.exp(2j * theta(t, r['Q'], r['mu'])) * v / v.conjugate())
    print(k, 'eps', r['eps'], flush=True, end=' '); print( 'ratios', [complex(round(z.real, 9), round(z.imag, 9)) for z in out])
