"""Turn results_ROUND.jsonl into the verdict table (markdown) and verdicts_ROUND.json.

For every statistic s and window:
  x        value on the true primes of the whole window
  T        theory value on the whole window (if a formula exists)
  z_T      block test: r_b = x_b - T_b on the K blocks, z_T = mean(r) / (sd(r)/sqrt(K))
  rel      (x - T) / |T|
  cramer   mean and sd over the coin-flip draws; z_C = (x - mean) / sd
  shuffle  the same for the shuffled gaps
Verdicts (fixed in preregister/round1.txt):
  structure  |z_C| > 3   (the primes differ from coins; for order statistics also |z_shuffle| > 3)
  theory     fits if |z_T| <= 3, misses otherwise
"""
import json, sys, os, math
from collections import defaultdict
import numpy as np

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')


def load(rnd):
    rows = [json.loads(l) for l in open(os.path.join(D, 'results_%s.jsonl' % rnd))]
    return rows


def analyse(rows, win):
    P = [r for r in rows if r['window'] == win and r['gen'] == 'primes']
    whole = [r for r in P if r['block'] == -1]
    if not whole:
        return {}
    whole = whole[-1]
    blocks = [r for r in P if r['block'] >= 0][-10:]
    nulls = defaultdict(list)
    for r in rows:
        if r['window'] == win and r['gen'] in ('cramer', 'shuffle'):
            nulls[r['gen']].append(r['stats'])
    out = {}
    for s, x in whole['stats'].items():
        d = {'x': x}
        T = whole['theory'].get(s)
        if T is not None:
            d['T'] = T
            d['rel'] = (x - T) / abs(T) if T else None
            if not s.startswith('A.race'):
                rb = np.array([b['stats'][s] - b['theory'][s] for b in blocks])
                se = rb.std(ddof=1) / math.sqrt(len(rb))
                d['zT'] = float(rb.mean() / se) if se > 0 else None
                d['r_block_mean'] = float(rb.mean()); d['r_block_se'] = float(se)
        for g, L in nulls.items():
            v = np.array([st[s] for st in L])
            d[g] = [float(v.mean()), float(v.std(ddof=1)), len(v)]
            d['z_' + g] = float((x - v.mean()) / v.std(ddof=1)) if v.std(ddof=1) > 0 else None
        zc = d.get('z_cramer')
        d['structure'] = None if zc is None else bool(abs(zc) > 3)
        if 'zT' in d and d['zT'] is not None:
            d['theory'] = 'fits' if abs(d['zT']) <= 3 else 'misses'
        out[s] = d
    return out


def fmt(v, n=4):
    if v is None:
        return ''
    if isinstance(v, float):
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e4):
            return '%.3e' % v
        return ('%.' + str(n) + 'f') % v
    return str(v)


def table(res):
    lines = ['| statistic | primes | theory | rel. diff | z (blocks) | coins | z coins | shuffled | z shuffled | verdict |',
             '|---|---|---|---|---|---|---|---|---|---|']
    for s, d in res.items():
        c = d.get('cramer'); sh = d.get('shuffle')
        lines.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
            s, fmt(d['x']), fmt(d.get('T')), fmt(d.get('rel'), 3) if d.get('rel') is not None else '',
            fmt(d.get('zT'), 1), (fmt(c[0]) + ' ± ' + fmt(c[1])) if c else '', fmt(d.get('z_cramer'), 1),
            (fmt(sh[0]) + ' ± ' + fmt(sh[1])) if sh else '', fmt(d.get('z_shuffle'), 1),
            d.get('theory', '') or ('structure' if d.get('structure') else '')))
    return '\n'.join(lines)


if __name__ == '__main__':
    rnd = sys.argv[1]
    rows = load(rnd)
    allres = {}
    for win in ('train', 'A', 'B'):
        res = analyse(rows, win)
        if res:
            allres[win] = res
            print('\n### window', win, '\n')
            print(table(res))
    json.dump(allres, open(os.path.join(D, 'verdicts_%s.json' % rnd), 'w'), indent=1)
