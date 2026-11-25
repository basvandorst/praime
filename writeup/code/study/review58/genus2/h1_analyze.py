"""H1: convergence of c = prod G_p(1/p) from own local data (cown/*.npz, exact p <= 1000, first order to 1e5).
Rules: raw partial product c(P); registered rule LC(P, nu): exp of the 1/p-weighted mean over primes P' in [P/10, P] of
log c(P') + nu/(sqrt P' ln P'); a smooth (Riesz) mean: exp of sum_p log G_p (1 - log p/log P)^2-free variant not used.
Out-of-sample: the rule at P = 1e4 / 5e4 against the raw product and the rule at 1e5. Slopes with old and new c."""
import sys, os, math, glob, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load, slope
rows = {r['key']: r for r in load()}
def LC(p, cum, P, nu):
    sel = (p >= P // 10) & (p <= P); w = 1.0 / p[sel]
    return float(np.sum((cum[sel] + nu / (np.sqrt(p[sel]) * np.log(p[sel]))) * w) / w.sum())
res = {}
for fn in sorted(glob.glob('cown/*.npz')):
    z = np.load(fn); key = os.path.basename(fn)[:-4]; p = z['p']; lg = z['logG']; sg = float(np.prod(z['sgn'])) if 'sgn' in z.files else 1.0
    cum = np.cumsum(lg)
    raw = {P: float(cum[p <= P][-1]) for P in (1000, 3000, 10000, 30000, 50000, 100000)}
    lc = {(P, nu): LC(p, cum, P, nu) for P in (10000, 50000, 100000) for nu in (0, 1, 2, 3)}
    # second-order term neglected beyond P0: exact minus first-order on 300 < p <= 1000 (good p), scaled by p^2
    e1, e2 = z['e1'], z['e2']; N = int(z['N'])
    m = (p > 300) & (p <= 1000) & np.array([N % int(q) != 0 for q in p])
    d2 = (lg[m] - np.log(1 - (e1[m] ** 2 - 1) / p[m])) * p[m] ** 2
    tail2 = float(np.sqrt(np.mean(d2 ** 2) * np.sum(1.0 / p[p > 1000].astype(float) ** 4)))   # random-walk size of the neglected sum
    sys2 = float(np.mean(d2) * np.sum(1.0 / p[p > 1000].astype(float) ** 2))                  # if its mean persisted
    # normalised tail fluctuation: (log c(1e5) - log c(P)) sqrt(P) ln P on P in [1e3, 3e4]
    sel = (p >= 1000) & (p <= 30000)
    K = (cum[-1] - cum[sel]) * np.sqrt(p[sel]) * np.log(p[sel])
    res[key] = dict(sign=sg, raw=raw, lc=lc, tail2=tail2, sys2=sys2, Krms=float(np.sqrt(np.mean(K ** 2))), Kmean=float(K.mean()),
                    M2=float(np.mean(e1[(p > 100)] ** 2)), N=N)
keys = sorted(res, key=lambda k: -abs(rows[k]['c']))
print('%-13s %5s %4s %8s | %8s %8s %8s | %8s %8s %8s | %7s %7s | %6s %6s' % ('key', 'N', 'eps', 'c_old', 'raw1e4', 'raw5e4', 'raw1e5',
      'LC5e4n2', 'LC1e5n2', 'LC1e5n0', 'new/old', 'oos', 'Krms', 'tail2'))
for k in keys:
    r = res[k]; s = r['sign']; co = rows[k]['c']
    new = s * math.exp(r['lc'][(100000, 2)])
    oos = s * math.exp(r['lc'][(50000, 2)]) / (s * math.exp(r['raw'][100000])) - 1      # registered rule at 5e4 vs raw at 1e5
    r['c_new'] = new
    print('%-13s %5d %+4.0f %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %+7.4f %+7.4f | %6.2f %.1e' % (k, r['N'], rows[k]['eps'], co,
          s * math.exp(r['raw'][10000]), s * math.exp(r['raw'][50000]), s * math.exp(r['raw'][100000]),
          s * math.exp(r['lc'][(50000, 2)]), new, s * math.exp(r['lc'][(100000, 0)]), new / co - 1, oos, r['Krms'], r['tail2']))
rel = np.array([res[k]['c_new'] / rows[k]['c'] - 1 for k in keys])
rel_rule = np.array([math.exp(res[k]['lc'][(50000, 2)]) / abs(rows[k]['c']) - 1 for k in keys])
oos = np.array([math.exp(res[k]['lc'][(50000, 2)] - res[k]['raw'][100000]) - 1 for k in keys])
oos2 = np.array([math.exp(res[k]['lc'][(50000, 2)] - res[k]['lc'][(100000, 2)]) - 1 for k in keys])
nus = np.array([[math.exp(res[k]['lc'][(100000, nu)] - res[k]['lc'][(100000, 2)]) - 1 for nu in (0, 1, 3)] for k in keys])
print('\nn = %d curves' % len(keys))
print('own rule at 5e4 (nu = 2) vs registered c: mean %+.5f, max |.| %.5f  (same rule, own local data)' % (rel_rule.mean(), np.abs(rel_rule).max()))
print('new c (rule at 1e5, nu = 2) / old c - 1: mean %+.5f, sd %.5f, min %+.5f, max %+.5f' % (rel.mean(), rel.std(), rel.min(), rel.max()))
print('out of sample: rule(5e4)/raw(1e5) - 1: rms %.5f, max |.| %.5f ; rule(5e4)/rule(1e5) - 1: rms %.5f, max %.5f' % (
      np.sqrt(np.mean(oos ** 2)), np.abs(oos).max(), np.sqrt(np.mean(oos2 ** 2)), np.abs(oos2).max()))
print('nu sensitivity at 1e5 (nu = 0, 1, 3 vs 2): max |.| %.5f %.5f %.5f' % tuple(np.abs(nus).max(0)))
Kr = np.array([res[k]['Krms'] for k in keys]); print('normalised tail K = (log c(1e5) - log c(P)) sqrt(P) ln P on [1e3, 3e4]: rms over curves %.2f, max %.2f'
      ' -> implied error of log c at P = 1e5: %.1e (rms), %.1e (max)' % (np.sqrt(np.mean(Kr ** 2)), Kr.max(), np.sqrt(np.mean(Kr ** 2)) / (math.sqrt(1e5) * math.log(1e5)), Kr.max() / (math.sqrt(1e5) * math.log(1e5))))
print('neglected second-order terms beyond p = 1000: random-walk size max %.1e, persistent-mean size max %.1e' % (max(res[k]['tail2'] for k in keys), max(abs(res[k]['sys2']) for k in keys)))
print('mean |lambda_p|^2 (100 < p <= 1e5): min %.4f max %.4f (USp(4): 1, sd of the mean ~ %.4f)' % (min(res[k]['M2'] for k in keys), max(res[k]['M2'] for k in keys), math.sqrt(2 / 9567)))
# correlations of the relative change with c, conductor, size of the tail correction (rule - raw at 5e4)
cc = np.array([rows[k]['c'] for k in keys]); NN = np.array([res[k]['N'] for k in keys]); tc = np.array([res[k]['lc'][(50000, 2)] - res[k]['raw'][50000] for k in keys])
for nm, x in (('c', cc), ('log N', np.log(NN)), ('tail correction', tc)):
    print('corr(new/old - 1, %s) = %+.3f' % (nm, np.corrcoef(rel, x)[0, 1]))
# slopes with old and new c (new c where computed)
R = list(rows.values())
def sl(S, which):
    return slope([r['G'] for r in S], [-r['eps'] * (res[r['key']]['c_new'] if (which == 'new' and r['key'] in res) else r['c']) * 150 / (2 * math.pi) for r in S])
for nm, S in (('all 100', R), ('eps+1', [r for r in R if r['eps'] > 0]), ('r151+152', [r for r in R if r['rnd'] != 'rep']), ('rep', [r for r in R if r['rnd'] == 'rep'])):
    print('slope %-9s old c %.4f +- %.4f   new c %.4f +- %.4f' % (nm, *sl(S, 'old'), *sl(S, 'new')))
json.dump({k: dict(c_old=rows[k]['c'], c_new=res[k]['c_new'], raw1e5=res[k]['sign'] * math.exp(res[k]['raw'][100000]),
                   rule5e4=res[k]['sign'] * math.exp(res[k]['lc'][(50000, 2)]), N=res[k]['N'], eps=rows[k]['eps'], M2=res[k]['M2']) for k in keys},
          open('c_new.json', 'w'), indent=1)
