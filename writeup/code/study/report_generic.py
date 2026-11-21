"""Block test of a theory file against a results file with lines {window, block, stats}.
usage: report_generic.py THEORY_JSON RESULTS_JSONL OUT_JSON"""
import json, sys, os, math
import numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
th = json.load(open(os.path.join(D, sys.argv[1])))
rows = [json.loads(l) for l in open(os.path.join(D, sys.argv[2]))]
out = {}
for win in th:
    P = [r for r in rows if r['window'] == win]
    whole = [r for r in P if r['block'] == -1][0]
    blocks = sorted([r for r in P if r['block'] >= 0], key=lambda r: r['block'])
    out[win] = {}
    print('\n### window %s\n' % win)
    print('| statistic | primes | theory | rel. diff | z (blocks) | verdict |\n|---|---|---|---|---|---|')
    for s, T in th[win]['-1'].items():
        if s not in whole['stats']:
            continue
        x = whole['stats'][s]
        rb = np.array([b['stats'][s] - th[win][str(b['block'])][s] for b in blocks])
        se = rb.std(ddof=1) / math.sqrt(len(rb))
        z = float(rb.mean() / se)
        o = {'x': x, 'T': T, 'rel': (x - T) / abs(T), 'zT': z, 'verdict': 'fits' if abs(z) <= 3 else 'misses'}
        for k in ('order1', 'order2'):
            if k in th[win]:
                o[k] = th[win][k][s]
        out[win][s] = o
        print('| %s | %.5g | %.5g | %+.2e | %.1f | %s |' % (s, x, T, o['rel'], z, o['verdict']))
json.dump(out, open(os.path.join(D, sys.argv[3]), 'w'), indent=1)
