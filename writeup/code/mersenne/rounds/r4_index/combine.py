"""Empirical p-values of every statistic of every derived sequence, and the family-wise correction
(Westfall-Young single-step min-p: the null distribution of the smallest p over all statistics, with leave-one-out
p-values of the null rows; row r of every W-derived sequence comes from the same synthetic list r).
Inputs: res_{d1,inv-d1,d2c,tt}.pkl (step rules), other.json, oeis.json.  Output: combined.txt
"""
import pickle, json, glob, os, numpy as np

D = np.load('data.npz')
NAMES = [str(x) for x in D['names']]
stats = []  # (seq, stat, real, null array)
# ---- step rules, all stages pooled ----
R = [pickle.load(open(f, 'rb')) for f in sorted(glob.glob('res_*.pkl'))]
stages = [r['stage'] for r in R]
blocks = R[0]['blocks']
Gtot = sum(r['G'] for r in R)
bestrules = {}
for nm, lo, hi in blocks:
    for k in ('exact', 'run0', 'runany'):
        v = np.max([r['best'][k][lo:hi] for r in R], axis=0)
        stats.append((nm, 'step_' + k, v[0], v[1:]))
    bi = np.max([r['oos']['best_in'][lo:hi] for r in R], axis=0)
    mx = np.max([np.where(r['oos']['best_in'][lo:hi] == bi, r['oos']['max'][lo:hi], -1) for r in R], axis=0)
    n = sum(np.where(r['oos']['best_in'][lo:hi] == bi, r['oos']['n'][lo:hi], 0) for r in R)
    sm = sum(np.where(r['oos']['best_in'][lo:hi] == bi, r['oos']['sum'][lo:hi], 0) for r in R)
    stats.append((nm, 'step_oosmax', mx[0], mx[1:]))
    stats.append((nm, 'step_oosmean', sm[0] / n[0], (sm / n)[1:]))
    # real best rules
    e = max(r['realbest'][nm][0] for r in R)
    bestrules[nm] = dict(exact=e, rules=[x for r in R if r['realbest'][nm][0] == e for x in r['realbest'][nm][1]][:6],
                         best_in=int(bi[0]), n_ties=int(n[0]), oos_max=int(mx[0]), oos_mean=float(sm[0] / n[0]),
                         oos_rules=[x for r in R if r['realoos'][nm][0] == bi[0] for x in r['realoos'][nm][1] if x[1] > 0][:5])
for fn in ('other.json', 'oeis.json'):
    if os.path.exists(fn):
        d = json.load(open(fn))
        for nm in NAMES:
            for k, rec in d[nm].items():
                stats.append((nm, k, float(rec['real']), np.array(rec['null'], float)))
nn = min(len(s[3]) for s in stats)
stats = [(a, b, float(c), np.asarray(d, float)[:nn]) for a, b, c, d in stats]


def p_real(r, v):
    return (1 + np.sum(v >= r)) / (len(v) + 1)


def p_loo(v):
    n = v.size; s = np.sort(v)
    ge = n - np.searchsorted(s, v, 'left') - 1
    return (1 + ge) / n


pr = np.array([p_real(s[2], s[3]) for s in stats])
PL = np.array([p_loo(s[3]) for s in stats]).T     # (nnull, nstat)
minp = PL.min(1)
adj = (1 + np.sum(minp[None, :] <= pr[:, None], axis=1)) / (nn + 1)
out = [f'{len(stats)} statistics over {len(NAMES)} derived sequences, {nn} null rows each; step grammar: '
       f'{Gtot:,} maps ({", ".join(stages)})',
       f'null distribution of the smallest p over all statistics: median {np.median(minp):.4f}, 5% {np.quantile(minp, .05):.4f}',
       f'real smallest raw p = {pr.min():.4f} ({stats[int(pr.argmin())][0]}/{stats[int(pr.argmin())][1]}), '
       f'family-wise adjusted p = {adj[int(pr.argmin())]:.3f}', '']
for nm in NAMES:
    idx = [i for i, s in enumerate(stats) if s[0] == nm]
    mp_seq = PL[:, idx].min(1)
    i0 = idx[int(np.argmin(pr[idx]))]
    adj_seq = (1 + np.sum(mp_seq <= pr[i0])) / (nn + 1)
    b = bestrules[nm]
    out.append(f'== {nm}: {len(idx)} stats; smallest raw p {pr[i0]:.4f} ({stats[i0][1]}: real {stats[i0][2]:g}, '
               f'null median {np.median(stats[i0][3]):g}); adjusted within sequence {adj_seq:.3f}, over all {adj[i0]:.3f}')
    out.append(f'   best step rule: {b["exact"]} exact transitions, e.g. ' +
               '; '.join(f'{r} at {h}' for r, h in b['rules'][:3]))
    out.append(f'   step OOS: best in-sample {b["best_in"]} hits ({b["n_ties"]} tied rules), max OOS {b["oos_max"]}, '
               f'mean OOS {b["oos_mean"]:.3f}' + (f'; rules with OOS hits: {b["oos_rules"]}' if b['oos_rules'] else ''))
    for i in idx:
        v = stats[i][3]
        out.append(f'   {stats[i][1]:14s} real {stats[i][2]:>9g}  null med {np.median(v):>9g} 95% {np.quantile(v, .95):>9g}  '
                   f'p {pr[i]:.4f}  fw {adj[i]:.3f}')
order = np.argsort(pr)
out.append('\nten smallest raw p over everything:')
for i in order[:10]:
    out.append(f'   {stats[i][0]:6s} {stats[i][1]:14s} real {stats[i][2]:g} null median {np.median(stats[i][3]):g} '
               f'p {pr[i]:.4f} family-wise {adj[i]:.3f}')
oo = [i for i, s in enumerate(stats) if 'oos' in s[1]]
out.append('\nout-of-sample statistics (real vs null):')
for i in oo:
    out.append(f'   {stats[i][0]:6s} {stats[i][1]:14s} real {stats[i][2]:g}  null mean {np.mean(stats[i][3]):.3f}  '
               f'share of nulls >= real {np.mean(stats[i][3] >= stats[i][2]):.3f}  p {pr[i]:.4f}')
open('combined.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
