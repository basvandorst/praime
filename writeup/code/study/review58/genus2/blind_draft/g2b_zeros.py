"""DRAFT. Zeros to T = 150 for every L-function of predictions_g2b.json with the registered routine of round 152
(bz4_lib.zeros, read-only import, completeness rule unchanged), the root number set from lfunrootres before lfuninit
(lfungenus2 leaves it unset; review 58, curve 1403). Resumable: one npz per L in zeros_g2b/.
usage: (ulimit -v 1500000; timeout 600 nice .venv/bin/python -I g2b_zeros.py BUDGET_SECONDS)   # repeat until ALL DONE"""
import sys, os, json, time
B4 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4'
sys.path.insert(0, B4); sys.path.insert(0, os.path.dirname(B4))
import numpy as np, cypari2, bz4_lib
HERE = os.path.dirname(os.path.abspath(__file__)); budget = float(sys.argv[1]); t0 = time.time()
_build = bz4_lib.build
def _build_eps(pari_, spec):
    if spec[0] != 'g2': return _build(pari_, spec)
    return pari_('my(L = lfungenus2([%s,%s])); L[6] = lfunrootres(L)[3]; L' % (spec[1], spec[2]))
bz4_lib.build = _build_eps
pari = cypari2.Pari(); pari.allocatemem(1100 * 10 ** 6)
P = json.load(open(os.path.join(HERE, 'predictions_g2b.json')))['L']; os.makedirs(os.path.join(HERE, 'zeros_g2b'), exist_ok=True)
left = 0
for key, e in sorted(P.items()):
    fn = os.path.join(HERE, 'zeros_g2b', key + '.npz')
    if os.path.exists(fn): continue
    if time.time() - t0 > budget - 200: left += 1; continue
    g, info = bz4_lib.zeros(pari, e['spec'], e['T2'])
    np.savez_compressed(fn, g=g, info=json.dumps(info)); print(key, info['status'], info['n'], '%.0f s' % (time.time() - t0), flush=True)
print('ALL DONE' if left == 0 else 'left %d' % left)
