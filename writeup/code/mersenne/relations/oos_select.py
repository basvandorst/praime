"""Out-of-sample check for the per-term statistics of angles 3 and 4.

Selection on terms 1..25: for every (base, feature) of angle 4, and every residue
class r mod 2^L (L = 2..4, odd r) of the terms' last bits, compute the empirical p of
the real first-25 value against the nulls' first-25 values.  Pick the smallest p.
Test: the same statistic on terms 26..52 against the nulls' terms 26..52.
Also reported: the full-data winner (run_23) split into its two halves.
"""
import json
import numpy as np
from angle4_bases import base_feats, FEATS, BASES
from common import emp_p

D = json.load(open('nulls.json'))
A = D['real']


def feats(seq):
    out = {}
    for b in BASES:
        f = base_feats(seq, b)
        for k in FEATS:
            out['%s_%d' % (k, b)] = f[k]
    for L in range(2, 5):
        for r in range(1, 1 << L, 2):
            out['res_%d_mod_%d' % (r, 1 << L)] = sum(1 for x in seq if x >= (1 << L) and x % (1 << L) == r)
    return out


def side(name):
    if name.startswith('res_'):
        return 'two'
    return FEATS[name.split('_')[0]]


for fam in ('a', 'b'):
    first = [feats(s[:25]) for s in D[fam]]
    second = [feats(s[25:]) for s in D[fam]]
    rf, rs = feats(A[:25]), feats(A[25:])
    ps = {k: emp_p(rf[k], [f[k] for f in first], side(k)) for k in rf}
    best = min(ps, key=ps.get)
    p_out = emp_p(rs[best], [f[best] for f in second], side(best))
    print('null (%s): chosen on terms 1..25: %s (real %.4g, p_in %.4f); terms 26..52: real %.4g, null median %.4g, p_out %.4f'
          % (fam, best, rf[best], ps[best], rs[best], np.median([f[best] for f in second]), p_out))
    for k in ('run_23', 'res_1_mod_4', 'res_1_mod_8', 'pal_2'):
        print('   %s: first half real %.4g p %.4f | second half real %.4g p %.4f'
              % (k, rf[k], emp_p(rf[k], [f[k] for f in first], side(k)),
                 rs[k], emp_p(rs[k], [f[k] for f in second], side(k))))
