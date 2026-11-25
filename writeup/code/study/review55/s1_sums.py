"""Review of round 151: own recomputation of every primitive sum from the zero files (own theta), completeness
diagnostics, pooled slopes with registered weights, jackknife, per-degree pools.  Reads only data files."""
import json, math, glob, os, sys
import numpy as np
from scipy.special import loggamma
import mpmath as mp

B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'
P = json.load(open(B + '/predictions_blind.json'))
RB = json.load(open(B + '/results_blind.json'))
OUT = os.path.dirname(os.path.abspath(__file__))

def theta(t, Q, mu):
    t = np.asarray(t, float)
    s = t / 2 * math.log(Q / math.pi ** len(mu))
    for m in mu:
        s = s + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return s

def theta_mp(t, Q, mu):
    t = mp.mpf(t)
    s = t / 2 * mp.log(mp.mpf(Q) / mp.pi ** len(mu))
    for m in mu:
        s += mp.im(mp.loggamma((mp.mpf('0.5') + m + 1j * t) / 2))
    return float(s)

def load(key):
    z = np.load(os.path.join(B, 'zeros_blind', key + '.npz'))
    return z['g'], json.loads(str(z['info']))

LINES = {'2': 2.0, '1/2': 0.5, '3': 3.0, '1/3': 1 / 3}
rows = {}
maxdiff = 0.0
for key, e in P['L'].items():
    g, info = load(key)
    T2 = e['T2']; Q = e['Q']; mu = e['mu']
    # theta spot check against mpmath
    for t in (g[0], g[len(g) // 2], g[-1]):
        maxdiff = max(maxdiff, abs(theta(t, Q, mu) - theta_mp(t, Q, mu)))
    gg = g[(g > 0) & (g <= T2)]
    th = theta(gg, Q, mu)
    ph = np.exp(2j * th)
    w = np.sin(math.pi * gg / T2) ** 2
    G = complex(ph.sum()); Gh = complex((w * ph).sum())
    Gl = {nm: complex((ph * np.exp(-1j * gg * math.log(r))).sum()) for nm, r in LINES.items()}
    # halves
    lo = gg <= T2 / 2
    Glo = complex(ph[lo].sum()); Ghi = complex(ph[~lo].sum())
    eps = complex(*e['eps']); c = e['c']
    pred_mine = -eps * c * T2 / (2 * math.pi)
    pred = complex(*e['pred'])
    # completeness diagnostics on all zeros in the file
    r = np.arange(1, len(g) + 1) - theta(g, Q, mu) / math.pi
    m = max(5, int(0.15 * len(g)))
    u = theta(g, Q, mu) / math.pi
    gaps = np.diff(u)
    rows[key] = dict(cls=e['cls'], deg=e['degree'], Q=Q, mu=mu, T2=T2, eps=[eps.real, eps.imag], c=c, M=e['M'],
                     sd_along=e['sd']['along'], sd_across=e['sd']['across'], sd_hann=e['sd']['along_hann'], sd_line=e['sd']['line_along'],
                     pred=[pred.real, pred.imag], pred_mine_diff=abs(pred - pred_mine), pred_hann=e['pred_hann'],
                     lines_pred={nm: e['lines'][nm]['pred'] for nm in LINES},
                     G=[G.real, G.imag], Gh=[Gh.real, Gh.imag], Gl={k: [v.real, v.imag] for k, v in Gl.items()},
                     Glo=[Glo.real, Glo.imag], Ghi=[Ghi.real, Ghi.imag], n=int(len(gg)), n_file=int(len(g)),
                     status=info['status'], m0=info['m0'], c_exp_reg=info['c_exp'],
                     r_mean=float(r.mean()), drift=float(r[-m:].mean() - r[:m].mean()),
                     r_top=float(r[-m:].mean()), max_gap=float(gaps.max()), n_gap_gt2=int((gaps > 2.0).sum()),
                     min_gap=float(gaps.min()), n_gap_lt01=int((gaps < 0.1).sum()), n_gap_lt02=int((gaps < 0.2).sum()),
                     n_gaps=int(len(gaps)), dup=int((np.diff(g) < 1e-6).sum()), g_first=float(g[0]), g_last=float(g[-1]))
    if key in RB['per_L']:
        Gr = complex(*RB['per_L'][key]['G'])
        rows[key]['G_repo_diff'] = abs(Gr - G)
print('theta: max |scipy - mpmath| =', maxdiff)
print('pred formula: max |pred - (-eps c T2/2pi)| =', max(r['pred_mine_diff'] for r in rows.values()))
print('max |G_mine - G_repo| =', max(r.get('G_repo_diff', 0) for r in rows.values()))
json.dump(rows, open(OUT + '/s1_rows.json', 'w'), indent=0)

# ------------- completeness summary per class
print('\ncompleteness (all zeros in file): r_mean - c_exp(reg), drift, max unfolded gap, gaps>2, gaps<0.1 (GUE expectation), dups')
def gue_cdf(x):  # Wigner surmise GUE: p(s) = 32/pi^2 s^2 exp(-4 s^2/pi)
    return float(mp.quad(lambda s: 32 / mp.pi ** 2 * s ** 2 * mp.exp(-4 * s ** 2 / mp.pi), [0, x]))
p01 = gue_cdf(0.1); p02 = gue_cdf(0.2); pgt2 = 1 - gue_cdf(2.0)
for cls in ['T1', 'T2', 'T3', 'T4a', 'T4b']:
    ks = [k for k in rows if rows[k]['cls'] == cls]
    dr = [rows[k]['r_mean'] - rows[k]['c_exp_reg'] for k in ks]
    ng = sum(rows[k]['n_gaps'] for k in ks)
    print('%-4s n=%2d  |r_mean-c_exp| max %.3f (excl. rejected %.3f)  |drift| max %.3f  max gap %.2f  #gap>2 %d (Wigner %.1f)  #gap<0.1 %d (Wigner %.1f)  #gap<0.2 %d (Wigner %.1f)  min gap %.4f  dups %d' % (
        cls, len(ks), max(abs(x) for x in dr), max(abs(rows[k]['r_mean'] - rows[k]['c_exp_reg']) for k in ks if rows[k]['status'] == 'COMPLETE'),
        max(abs(rows[k]['drift']) for k in ks), max(rows[k]['max_gap'] for k in ks), sum(rows[k]['n_gap_gt2'] for k in ks), pgt2 * ng,
        sum(rows[k]['n_gap_lt01'] for k in ks), p01 * ng, sum(rows[k]['n_gap_lt02'] for k in ks), p02 * ng,
        min(rows[k]['min_gap'] for k in ks), sum(rows[k]['dup'] for k in ks)))
for k in rows:
    if abs(rows[k]['r_mean'] - rows[k]['c_exp_reg']) > 0.05:
        print('   offset:', k, rows[k]['status'], 'r_mean %.4f c_exp %.4f diff %.4f' % (rows[k]['r_mean'], rows[k]['c_exp_reg'], rows[k]['r_mean'] - rows[k]['c_exp_reg']))
