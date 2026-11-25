"""Pooled slopes with the registered weights, delete-one jackknife over L-functions, per family and per degree; noise model
checked by chi^2 of the per-L 'along' residuals; amplitude regression (beta) per subset; c = 1 and alternatives."""
import json, math, os, sys
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(D + '/s1_rows.json'))
# the two chi mod 17 with their (complete) zero sets
for k in ('T1_17_5', 'T1_17_7'):
    s = json.load(open(D + '/s2_%s.json' % k)); R[k]['G'] = s['G']; R[k]['Gh'] = s['Gh']; R[k]['status'] = 'COMPLETE(review)'
def cx(v): return complex(v[0], v[1])

def pooled(keys, which='r1'):
    num = 0.0; numc = 0j; var = 0.0
    for k in keys:
        r = R[k]
        if which == 'r1': pr, G, sa = cx(r['pred']), cx(r['G']), r['sd_along']
        elif which == 'hann': pr, G, sa = cx(r['pred_hann']), cx(r['Gh']), r['sd_hann']
        elif which == 'lo': pr, G, sa = cx(r['pred']) / 2, cx(r['Glo']), r['sd_along'] / math.sqrt(2)
        elif which == 'hi': pr, G, sa = cx(r['pred']) / 2, cx(r['Ghi']), r['sd_along'] / math.sqrt(2)
        num += abs(pr) ** 2; numc += G * pr.conjugate(); var += abs(pr) ** 2 * sa ** 2
    return numc.real / num, math.sqrt(var) / num, numc.imag / num

def jack(keys, which='r1'):
    n = len(keys); full = pooled(keys, which)[0]
    loo = np.array([pooled([k for k in keys if k != kk], which)[0] for kk in keys])
    se = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    return se, loo

def chi2_along(keys, b):
    # residual of Re(G conj pr)/|pr| - b|pr| in units of sd_along
    z = []
    for k in keys:
        r = R[k]; pr = cx(r['pred']); G = cx(r['G'])
        z.append(((G * pr.conjugate()).real / abs(pr) - b * abs(pr)) / r['sd_along'])
    z = np.array(z); return float(np.sum(z ** 2)), len(z), z

reg = [k for k in R if R[k]['status'] == 'COMPLETE']
all40 = list(R)
fam = {c: [k for k in reg if R[k]['cls'] == c] for c in ['T1', 'T2', 'T3', 'T4a', 'T4b']}
def line(name, keys, which='r1'):
    b, se, im = pooled(keys, which); js, loo = jack(keys, which)
    c2, n, z = chi2_along(keys, b) if which == 'r1' else (float('nan'), len(keys), None)
    print('%-34s n=%3d  b = %.4f  model se %.4f  jackknife se %.4f  (b-1)/se_model %+.2f  (b-1)/se_jack %+.2f  chi2/n along %.2f  min/max LOO %.4f/%.4f' % (
        name, len(keys), b, se, js, (b - 1) / se, (b - 1) / js, c2 / n, loo.min(), loo.max()))
    return b, se, js
print('=== registered set (128) ===')
line('ALL (D2)', reg)
for c, ks in fam.items(): line('family ' + c, ks)
print('\n=== with the two chi mod 17 (complete; rejected by the arg-branch bug) ===')
line('ALL 130', all40)
line('T1 40', [k for k in all40 if R[k]['cls'] == 'T1'])
print('\n=== per degree ===')
for name, ks in [('degree 1 (T1)', fam['T1']), ('degree 2 (T2)', fam['T2']), ('degree 3 (T3)', fam['T3']),
                 ('degree 4 (T4a+T4b)', fam['T4a'] + fam['T4b']), ('degree >= 2', fam['T2'] + fam['T3'] + fam['T4a'] + fam['T4b']),
                 ('degree >= 3', fam['T3'] + fam['T4a'] + fam['T4b'])]:
    line(name, ks)
print('\n=== Hann ===')
for name, ks in [('ALL', reg), ('degree >= 2', fam['T2'] + fam['T3'] + fam['T4a'] + fam['T4b']), ('degree >= 3', fam['T3'] + fam['T4a'] + fam['T4b']),
                 ('T4a', fam['T4a']), ('T4b', fam['T4b']), ('degree 4', fam['T4a'] + fam['T4b'])]:
    line('Hann ' + name, ks, 'hann')
print('\n=== lower half / upper half of each window (sharp) ===')
for name, ks in [('T1', fam['T1']), ('T2', fam['T2']), ('T3', fam['T3']), ('T4a', fam['T4a']), ('T4b', fam['T4b']), ('degree 4', fam['T4a'] + fam['T4b'])]:
    blo, slo, _ = pooled(ks, 'lo'); bhi, shi, _ = pooled(ks, 'hi')
    print('%-10s lower half %.3f ± %.3f   upper half %.3f ± %.3f' % (name, blo, slo, bhi, shi))

# ---- per-degree pooled excluding objects with tiny |c| leverage? (none: weights |pred|^2 already)
# ---- weight shares
tot = sum(abs(cx(R[k]['pred'])) ** 2 for k in reg)
for c, ks in fam.items():
    print('weight share %-4s %.3f' % (c, sum(abs(cx(R[k]['pred'])) ** 2 for k in ks) / tot))

# ---- D5 regression of c_hat on [1, c-1] per subset; and chi^2 of alternatives
def d5(keys):
    X = []; Y = []; W = []
    for k in keys:
        r = R[k]; eps = cx(r['eps']); unit = r['T2'] / (2 * math.pi)
        y = (cx(r['G']) * (-eps).conjugate()).real / abs(eps) / unit; s = r['sd_along'] / unit
        X.append([1.0, r['c'] - 1.0]); Y.append(y); W.append(1 / s ** 2)
    X = np.array(X); Y = np.array(Y); W = np.array(W)
    A = X.T @ (W[:, None] * X); bv = np.linalg.solve(A, X.T @ (W * Y)); cov = np.linalg.inv(A)
    # proportional fit y = beta c (no intercept)
    cc = X[:, 1] + 1; bp = np.sum(W * cc * Y) / np.sum(W * cc * cc); sbp = 1 / math.sqrt(np.sum(W * cc * cc))
    # jackknife for beta
    n = len(keys); loo = []
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        Ai = X[m].T @ (W[m][:, None] * X[m]); loo.append(np.linalg.solve(Ai, X[m].T @ (W[m] * Y[m]))[1])
    loo = np.array(loo); js = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    chi = lambda pred: float(np.sum(W * (Y - pred) ** 2))
    return dict(alpha=bv[0], beta=bv[1], se_a=math.sqrt(cov[0, 0]), se_b=math.sqrt(cov[1, 1]), se_b_jack=js, prop=bp, se_prop=sbp,
                chi_law=chi(cc), chi_c1=chi(np.ones(n)), chi_0=chi(np.zeros(n)), n=n, Y=Y, c=cc, W=W)
print('\n=== amplitude regression c_hat = alpha + beta (c - 1) ===')
for name, ks in [('ALL', reg), ('degree >= 2', fam['T2'] + fam['T3'] + fam['T4a'] + fam['T4b']), ('degree >= 3', fam['T3'] + fam['T4a'] + fam['T4b']),
                 ('T2', fam['T2']), ('T3', fam['T3']), ('degree 4', fam['T4a'] + fam['T4b'])]:
    d = d5(ks)
    print('%-12s n=%3d alpha %.3f ± %.3f  beta %.3f ± %.3f (jack %.3f)  proportional %.3f ± %.3f  chi2: law %.1f  c=1 %.1f  H0 %.1f' % (
        name, d['n'], d['alpha'], d['se_a'], d['beta'], d['se_b'], d['se_b_jack'], d['prop'], d['se_prop'], d['chi_law'], d['chi_c1'], d['chi_0']))
json.dump({k: dict(G=R[k]['G'], pred=R[k]['pred'], c=R[k]['c'], cls=R[k]['cls']) for k in R}, open(D + '/s3_dump.json', 'w'))
