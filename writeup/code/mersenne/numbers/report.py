"""Summarise res_*.pkl: real scores, null distributions, empirical p-values, best rules, OOS, MDL."""
import pickle, math, glob, numpy as np, sys
sys.path.insert(0, '.')
from seqs import REAL

KEYS = ('exact', 'run0', 'runany', 'pm1', 'pm2', 'rel1')
R = {}
for f in sorted(glob.glob('res_*.pkl')):
    d = pickle.load(open(f, 'rb')); R[d['stage']] = d
out = []
def pr(s=''):
    out.append(s); print(s)

lam = 51 / math.log(REAL[-1] / REAL[0])
def bits_next(i, inv=False, tt=False):  # cost of coding the target of transition i under the Poisson-in-ln null
    if tt:
        x, y = REAL[i + 1], REAL[i + 2]
    elif inv:
        x, y = REAL[i + 1], REAL[i]
    else:
        x, y = REAL[i], REAL[i + 1]
    lo, hi = min(x, y), max(x, y)
    return -math.log2(min(1.0, lam * (lo / hi) ** lam / y))


def group(stages):
    ds = [R[s] for s in stages if s in R]
    if not ds: return None
    best = {k: np.max([d['best'][k] for d in ds], axis=0) for k in KEYS}
    G = sum(d['G'] for d in ds)
    # OOS: rule(s) with the best in-sample exact count over the union of stages
    bi = np.max([d['oos']['best_in'] for d in ds], axis=0)
    n = sum(np.where(d['oos']['best_in'] == bi, d['oos']['n'], 0) for d in ds)
    sm = sum(np.where(d['oos']['best_in'] == bi, d['oos']['sum'], 0) for d in ds)
    mx = np.max([np.where(d['oos']['best_in'] == bi, d['oos']['max'], -1) for d in ds], axis=0)
    return dict(best=best, G=G, bi=bi, n=n, mean=sm / n, mx=mx, nA=ds[0]['nA'], T=ds[0]['T'], NIN=ds[0]['NIN'],
                cands=[c for d in ds for c in d['cands']])


def null_table(g, title):
    nA = g['nA']
    pr(f'\n=== {title}: {g["G"]:,} maps, T = {g["T"]} transitions')
    pr(f'{"metric":8s} {"real":>5s} | {"nullA med":>9s} {"95%":>4s} {"max":>4s} {"p":>6s} | {"nullB med":>9s} {"95%":>4s} {"max":>4s} {"p":>6s}')
    for k in KEYS:
        v = g['best'][k]; r = v[0]; line = f'{k:8s} {r:5d} |'
        for sl in (slice(1, 1 + nA), slice(1 + nA, None)):
            x = v[sl]; p = (1 + (x >= r).sum()) / (1 + len(x))
            line += f' {np.median(x):9.1f} {np.quantile(x, .95):4.0f} {x.max():4d} {p:6.3f} |'
        pr(line)
    # distribution of the exact maximum
    for nm, sl in (('A', slice(1, 1 + nA)), ('B', slice(1 + nA, None))):
        x = g['best']['exact'][sl]; vals, cnt = np.unique(x, return_counts=True)
        pr(f'  null{nm} best-exact histogram: ' + ', '.join(f'{a}:{b}' for a, b in zip(vals, cnt)))
    pr(f'  OOS (choose on transitions 1..{g["NIN"]}, count hits on the remaining {g["T"] - g["NIN"]}):')
    for nm, sl in (('real', slice(0, 1)), ('nullA', slice(1, 1 + nA)), ('nullB', slice(1 + nA, None))):
        pr(f'    {nm:5s}: best in-sample exact = {np.median(g["bi"][sl]):.1f} (median), ties = {np.median(g["n"][sl]):.0f}, '
           f'mean OOS hits of a tied rule = {np.mean(g["mean"][sl]):.3f}, max OOS hits over ties = {np.mean(g["mx"][sl]):.2f} (mean)'
           + (f', 95% {np.quantile(g["mx"][sl], .95):.0f}' if nm != 'real' else ''))


def top_rules(g, title, inv=False, tt=False, k=8):
    C = g['cands']; G = g['G']; T = g['T']
    def mdl(c):
        hits = c[2]
        gain = sum(bits_next(i + (1 if tt else 0) - (1 if tt else 0), inv, tt) for i in hits)
        cost = math.log2(G) + math.log2(math.comb(T, len(hits)))
        return gain - cost
    pr(f'\n--- best rules on the real sequence, {title} (r = exact, run0, runany, pm1, pm2, rel1, in, oos); MDL net bits = gain - log2 G - log2 C(T,k)')
    for key, j in (('exact', 0), ('run0', 1), ('runany', 2), ('pm1', 3), ('rel1', 5), ('oos', 7)):
        srt = sorted(C, key=lambda c: (c[1][j], c[1][0], -len(c[0])), reverse=True)[:k]
        if not srt: continue
        pr(f'  by {key}:')
        for c in srt[:k]:
            tg = [f'{REAL[i + (1 if tt else 0)] if not inv else REAL[i + 1]}->{REAL[i + (2 if tt else 1)] if not inv else REAL[i]}' for i in c[2]]
            pr(f'    {c[1]}  MDL {mdl(c):6.1f}  {c[0]}   hits: {" ".join(tg)}')
    best_mdl = max(C, key=mdl) if C else None
    if best_mdl:
        pr(f'  highest MDL net: {mdl(best_mdl):.1f} bits: {best_mdl[0]} {best_mdl[1]}')


fw = group(['d1', 'd2', 'd3'])
if fw:
    for st in ('d1', 'd2', 'd3'):
        if st in R: null_table(group([st]), f'forward stage {st}')
    null_table(fw, 'forward successor maps, all depths')
    top_rules(fw, 'forward')
if 'tt' in R:
    g = group(['tt']); null_table(g, 'two-term rules a_{n+1} = h(a_n, a_{n-1})'); top_rules(g, 'two-term', tt=True)
iv = group(['inv-d1', 'inv-d2', 'inv-d3'])
if iv:
    null_table(iv, 'inverse maps a_n = f(a_{n+1}), all depths'); top_rules(iv, 'inverse', inv=True)
open('report.txt', 'w').write('\n'.join(out) + '\n')
