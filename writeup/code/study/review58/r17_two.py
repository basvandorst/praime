"""Residual profile (cumulative Re(G - pred)) for curves with 2 | N, with a non-squarefree conductor, and the rest."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
A = {}
for line in open('r4_arith.jsonl'):
    o = json.loads(line); A[o['key']] = o
rows = [r for r in load() if r['cls'] == 'g2']
Ts = [5, 10, 20, 30, 45, 60, 75, 90, 105, 120, 135, 150]
def prof(S):
    D = np.array([[(np.exp(2j * theta(r['g'][(r['g'] > 0) & (r['g'] <= T)], r['Q'], r['mu'])).sum() + r['eps'] * r['c'] * T / (2 * math.pi)).real for T in Ts] for r in S])
    return D
for lab, f in (('2 | N', lambda b: any(p == 2 for p, e in b)), ('odd N', lambda b: not any(p == 2 for p, e in b)),
               ('exp>=2 somewhere', lambda b: any(e >= 2 for p, e in b)), ('squarefree odd', lambda b: all(e == 1 for p, e in b) and not any(p == 2 for p, e in b))):
    S = [r for r in rows if f(A[r['key']]['bad'])]
    D = prof(S)
    print('%-18s n=%2d ' % (lab, len(S)) + ' '.join('T%d:%+.1f(%.1f)' % (T, *jk_mean(D[:, j])) for j, T in enumerate(Ts)))
for r in rows:
    if any(p == 2 for p, e in A[r['key']]['bad']):
        print('  ', r['key'], A[r['key']]['bad'], A[r['key']]['loc'].get('2'), 'eps', r['eps'], 'c %.3f' % r['c'])
