"""Block test of a theory file against the primes of round 1 (the statistics of the primes do not change).
usage: report2.py THEORY_JSON [ROUND_RESULTS=r1]  -> markdown tables and verdicts_<theory>.json"""
import json, sys, os, math
import numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
th = json.load(open(os.path.join(D, sys.argv[1])))
rows = [json.loads(l) for l in open(os.path.join(D, 'results_%s.jsonl' % (sys.argv[2] if len(sys.argv) > 2 else 'r1')))]
out = {}
for win in th:
    P = [r for r in rows if r['window'] == win and r['gen'] == 'primes']
    whole = [r for r in P if r['block'] == -1][-1]
    blocks = sorted([r for r in P if r['block'] >= 0], key=lambda r: r['block'])[-10:]
    print('\n### window %s\n' % win)
    print('| statistic | primes | theory | rel. diff | z (blocks) | verdict |\n|---|---|---|---|---|---|')
    out[win] = {}
    for s, T in th[win]['-1'].items():
        if s not in whole['stats']:
            continue
        x = whole['stats'][s]
        rb = np.array([b['stats'][s] - th[win][str(b['block'])][s] for b in blocks])
        se = rb.std(ddof=1) / math.sqrt(len(rb))
        z = rb.mean() / se
        v = 'fits' if abs(z) <= 3 else 'misses'
        out[win][s] = {'x': x, 'T': T, 'rel': (x - T) / abs(T), 'zT': float(z), 'verdict': v}
        print('| %s | %.6g | %.6g | %+.2e | %.1f | %s |' % (s, x, T, (x - T) / abs(T), z, v))
json.dump(out, open(os.path.join(D, 'verdicts_' + sys.argv[1]), 'w'), indent=1)
