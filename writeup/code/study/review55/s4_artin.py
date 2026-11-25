"""Artin factors of the 30 cubic fields: rho-zeros = zeros of zeta_K minus the zeros of zeta (own: mpmath zetazero),
own theta, sharp/Hann/smoothed-edge sums, completeness of the rho-zero lists, jackknife over fields, noise."""
import json, math, os
import numpy as np, mpmath as mp
from scipy.special import loggamma
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'
D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json'))
def theta(t, Q, mu):
    t = np.asarray(t, float); s = t / 2 * math.log(Q / math.pi ** len(mu))
    for m in mu: s = s + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return s
zf = D + '/zeta_zeros_310.npy'
if not os.path.exists(zf):
    mp.mp.dps = 20; zz = []
    n = 1
    while True:
        g = float(mp.zetazero(n).imag)
        if g > 310: break
        zz.append(g); n += 1
    np.save(zf, np.array(zz))
zz = np.load(zf)
print('zeta zeros <= 300 (mpmath):', int((zz <= 300).sum()))
from scipy.special import erf
def Hsm(t, T2, Dl):  # indicator (0,T2] convolved with Gaussian of width Dl
    return 0.5 * (erf(t / Dl) - erf((t - T2) / Dl))
rows = {}
for key, e in P['products'].items():
    if e['cls'] != 'T5': continue
    z = np.load(B + '/zeros_blind/' + key + '.npz'); g = z['g']; info = json.loads(str(z['info']))
    T2 = e['T2']; r = e['rho']
    # match zeta zeros
    idx = np.searchsorted(g, zz)
    d = np.full(len(zz), np.inf)
    for j, (i, x) in enumerate(zip(idx, zz)):
        for ii in (i - 1, i):
            if 0 <= ii < len(g): d[j] = min(d[j], abs(g[ii] - x))
    zmatch = zz[(zz <= g[-1] - 0.5)]
    nm = int((d[: len(zmatch)] < 1e-6).sum())
    isz = np.zeros(len(g), bool)
    for x in zz:
        i = np.argmin(np.abs(g - x))
        if abs(g[i] - x) < 1e-6: isz[i] = True
    gr = g[~isz]
    # completeness of rho list: n - theta_rho/pi, expected 1/2 (eps = +1, no central zero assumed) ; K list: expected 1.5
    ur = np.arange(1, len(gr) + 1) - theta(gr, r['Q'], r['mu']) / math.pi
    uk = np.arange(1, len(g) + 1) - theta(g, e['Q'], e['mu']) / math.pi
    m = max(5, int(0.15 * len(gr)))
    grr = gr[(gr > 0) & (gr <= T2)]
    ph = np.exp(2j * theta(grr, r['Q'], r['mu']))
    w = np.sin(math.pi * grr / T2) ** 2
    R = complex(ph.sum()); Rh = complex((w * ph).sum())
    Rs = {Dl: complex((Hsm(grr, T2, Dl) * ph).sum()) for Dl in (5.0, 20.0)}
    lines = {nm_: complex((ph * np.exp(-1j * grr * math.log(rr))).sum()) for nm_, rr in (('2', 2.), ('1/2', .5), ('3', 3.), ('1/3', 1 / 3))}
    lo = grr <= T2 / 2
    pr = e['readings']['R_artin']
    rows[key] = dict(d=e['d'], Q=r['Q'], mu=r['mu'], c=r['c'], c_direct=r.get('c_direct'), pred=pr['pred'], pred_hann=pr['pred_hann'], sd=pr['sd'],
                     sd_hann=pr['sd_hann'], M=r['M'], R=[R.real, R.imag], Rh=[Rh.real, Rh.imag],
                     Rs={str(k): [v.real, v.imag] for k, v in Rs.items()}, Rlo=[complex(ph[lo].sum()).real], Rhi=[complex(ph[~lo].sum()).real],
                     lines={k: [v.real, v.imag] for k, v in lines.items()}, lines_pred={k: v['pred'][0] for k, v in pr['lines'].items()},
                     n_rho=int(len(grr)), n_K=int(len(g)), zeta_matched=nm, zeta_total=int(len(zmatch)),
                     rmean_rho=float(ur.mean()), drift_rho=float(ur[-m:].mean() - ur[:m].mean()), rmean_K=float(uk.mean()),
                     min_gap_rho=float(np.min(np.diff(theta(gr, r['Q'], r['mu']) / math.pi))), m0=info['m0'])
json.dump(rows, open(D + '/s4_rows.json', 'w'), indent=0)
ks = sorted(rows, key=lambda k: -rows[k]['d'])
print('key        d     c_rho   pred    R(sharp)  slope  Hann-slope  n_rho  zeta-match  rmean_rho drift  rmean_K')
for k in ks:
    r = rows[k]
    print('%-10s %5d  %6.3f %7.1f %8.1f  %6.3f  %6.3f   %4d  %d/%d  %+.3f %+.3f  %+.3f' % (k, r['d'], r['c'], r['pred'], r['R'][0],
          r['R'][0] / r['pred'] if r['pred'] else float('nan'), r['Rh'][0] / r['pred_hann'] if r['pred_hann'] else float('nan'),
          r['n_rho'], r['zeta_matched'], r['zeta_total'], r['rmean_rho'], r['drift_rho'], r['rmean_K']))
def slope(keys, meas='R', pk='pred', sk='sd'):
    num = sum(rows[k][meas][0] * rows[k][pk] for k in keys); den = sum(rows[k][pk] ** 2 for k in keys)
    var = sum(rows[k][pk] ** 2 * rows[k][sk] ** 2 for k in keys)
    return num / den, math.sqrt(var) / den
def jack(keys, **kw):
    loo = np.array([slope([k for k in keys if k != kk], **kw)[0] for kk in keys]); n = len(keys)
    return math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2)), loo
ks30 = list(rows); ks29 = [k for k in ks30 if rows[k]['c'] != 0]
for name, kw in (('sharp', dict(meas='R', pk='pred', sk='sd')), ('Hann', dict(meas='Rh', pk='pred_hann', sk='sd_hann'))):
    b, se = slope(ks30, **kw); js, loo = jack(ks30, **kw)
    print('Artin %-5s slope %.4f  model se %.4f  jackknife se %.4f  -> (b-1)/se_model %.2f, (b-1)/se_jack %.2f; LOO range %.4f..%.4f' % (
        name, b, se, js, (b - 1) / se, (b - 1) / js, loo.min(), loo.max()))
for Dl in ('5.0', '20.0'):
    num = sum(rows[k]['Rs'][Dl][0] * rows[k]['pred'] for k in ks30); den = sum(rows[k]['pred'] ** 2 for k in ks30)
    print('Artin smoothed-edge window (Gaussian %s) slope %.4f' % (Dl, num / den))
num = sum(rows[k]['Rlo'][0] * rows[k]['pred'] / 2 for k in ks30); den = sum((rows[k]['pred'] / 2) ** 2 for k in ks30)
num2 = sum(rows[k]['Rhi'][0] * rows[k]['pred'] / 2 for k in ks30)
print('Artin lower half (0,150] slope %.4f   upper half (150,300] slope %.4f' % (num / den, num2 / den))
# lines
for nm_ in ('2', '1/2', '3', '1/3'):
    num = sum(rows[k]['lines'][nm_][0] * rows[k]['lines_pred'][nm_] for k in ks30); den = sum(rows[k]['lines_pred'][nm_] ** 2 for k in ks30)
    print('Artin line %-3s slope %.4f' % (nm_, num / den))
# noise: along variance / M of the rho zeros
b, _ = slope(ks30)
va = np.mean([(rows[k]['R'][0] - b * rows[k]['pred']) ** 2 / rows[k]['M'] for k in ks30])
vc = np.mean([rows[k]['R'][1] ** 2 / rows[k]['M'] for k in ks30])
print('Artin Var/M along %.3f  across %.3f  (registered model along 0.077 -> sd %s)' % (va, vc, 'sqrt(0.077 M)'))
print('chi2/n along with registered sd: %.2f' % np.mean([((rows[k]['R'][0] - b * rows[k]['pred']) / rows[k]['sd']) ** 2 for k in ks30]))
# split by conductor prime / not; by c
isprime = lambda n: all(n % p for p in range(2, int(n ** .5) + 1))
for name, sub in (('|d| prime', [k for k in ks30 if isprime(-rows[k]['d'])]), ('|d| composite', [k for k in ks30 if not isprime(-rows[k]['d'])]),
                  ('c_rho < 1.5', [k for k in ks30 if rows[k]['c'] < 1.5]), ('c_rho >= 1.5', [k for k in ks30 if rows[k]['c'] >= 1.5])):
    b, se = slope(sub); js, _ = jack(sub)
    print('Artin %-14s n=%2d slope %.4f ± %.4f (jack %.4f)' % (name, len(sub), b, se, js))
