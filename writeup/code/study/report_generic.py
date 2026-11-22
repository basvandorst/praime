"""Block test of a theory file against a results file with lines {window, block, stats}.
usage: report_generic.py THEORY_JSON RESULTS_JSONL OUT_JSON

z (blocks) is a Student t with 9 degrees of freedom (the sd comes from 10 blocks): p_t9 is its two-sided p-value.
For shares the binomial z = (x - T)/sqrt(T(1-T)/N) on the whole window is given as well (N = number of trials:
gaps for W1, W2, L.*; gaps <= 300 for W3). Verdict as pre-registered: fits when |z (blocks)| <= 3."""
import json, sys, os, math
import numpy as np
from scipy.stats import t as student_t
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
th = json.load(open(os.path.join(D, sys.argv[1])))
rows = [json.loads(l) for l in open(os.path.join(D, sys.argv[2]))]


def trials(s, st):
    if 'n_gaps' in st and (s.startswith('W1') or s.startswith('W2')):
        return st['n_gaps']
    if 'n_gaps' in st and s.startswith('W3'):
        return st['n_gaps'] * st['W2.le300']
    if 'count' in st and s.startswith('L.'):
        return st['count'] - 1
    return None


out = {}
nt = nmiss = 0
for win in th:
    P = [r for r in rows if r['window'] == win]
    whole = [r for r in P if r['block'] == -1][0]
    blocks = sorted([r for r in P if r['block'] >= 0], key=lambda r: r['block'])
    out[win] = {}
    print('\n### window %s\n' % win)
    print('| statistic | primes | theory | rel. diff | z (blocks) | p (t9) | z (binomial) | verdict |')
    print('|---|---|---|---|---|---|---|---|')
    for s, T in th[win]['-1'].items():
        if s not in whole['stats']:
            continue
        x = whole['stats'][s]
        rb = np.array([b['stats'][s] - th[win][str(b['block'])][s] for b in blocks])
        se = rb.std(ddof=1) / math.sqrt(len(rb))
        z = float(rb.mean() / se)
        p = float(2 * student_t.sf(abs(z), len(rb) - 1))
        N = trials(s, whole['stats'])
        zb = float((x - T) / math.sqrt(T * (1 - T) / N)) if N else None
        o = {'x': x, 'T': T, 'rel': (x - T) / abs(T), 'zT': z, 'p_t9': p, 'z_binom': zb,
             'verdict': 'fits' if abs(z) <= 3 else 'misses'}
        nt += 1; nmiss += abs(z) > 3
        for k in ('order1', 'order2'):
            if k in th[win]:
                o[k] = th[win][k][s]
        out[win][s] = o
        print('| %s | %.5g | %.5g | %+.2e | %.1f | %.3f | %s | %s |' % (s, x, T, o['rel'], z, p,
              '%.1f' % zb if zb is not None else '', o['verdict']))
p3 = float(2 * student_t.sf(3, 9))
print('\n%d tests, %d with |t9| > 3; expected by chance %.1f (p = %.4f each)' % (nt, nmiss, nt * p3, p3))
out['_summary'] = {'tests': nt, 'misses': int(nmiss), 'expected_by_chance': nt * p3}
json.dump(out, open(os.path.join(D, sys.argv[3]), 'w'), indent=1)
