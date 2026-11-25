# Block jackknife over several block counts, control-frequency noise, half splits, pooled slopes.
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, numpy as np
from r1_sums import theta, own_data, gauss_eps, D, OUT
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
r1 = json.load(open(OUT + 'r1_sums.json'))
L = own_data()
for key in ('X5', 'X5d', 'X7', 'X7d'):
    q = 5 if key.startswith('X5') else 7
    L[key] = dict(Q=float(q), mu=pred[key]['mu'], eps=complex(*pred[key]['eps']), fam='GL1 complex character')

from fractions import Fraction
rats = set()
for a in range(1, 7):
    for b in range(1, 7):
        rats.add(Fraction(a, b))
lnr = np.array(sorted(np.log(float(f)) for f in rats))

rng = np.random.default_rng(1)
out = {}
rows = []
for key, v in L.items():
    P = pred[key]
    T2 = float(P['T2'])
    g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']
    g = np.sort(g[(g > 0) & (g <= T2)])
    th = theta(g, v['Q'], v['mu'])
    e = np.exp(2j * th)
    p = -v['eps'] * P['c'] * T2 / (2 * np.pi)
    ph = p / abs(p)
    G = e.sum()
    b = (G / p).real
    # blocks
    jk = {}
    for K in (3, 4, 6, 8, 12, 20, 30):
        edges = np.linspace(0, T2, K + 1)
        idx = np.searchsorted(edges, g, side='left') - 1
        Gk = np.array([e[idx == k].sum() for k in range(K)])
        bk = (Gk / (p / K)).real
        se = bk.std(ddof=1) / np.sqrt(K)
        # delete-one jackknife of the total slope
        loo = np.array([((G - Gk[k]) / (p * (K - 1) / K)).real for k in range(K)])
        sej = np.sqrt((K - 1) / K * ((loo - loo.mean()) ** 2).sum())
        jk[K] = (float(se), float(sej))
    # halves
    lo = e[g <= T2 / 2].sum(); hi = e[g > T2 / 2].sum()
    b_lo = (lo / (p / 2)).real; b_hi = (hi / (p / 2)).real
    # control frequencies (sharp window and Hann window), avoiding all rationals a/b <= 30 by 3*2pi/T2
    w = np.sin(np.pi * g / T2) ** 2
    us = []
    while len(us) < 400:
        u = rng.uniform(0.15, 3.0) * rng.choice([-1, 1])
        if np.min(np.abs(lnr - u)) > 1.5 * 2 * np.pi / T2:
            us.append(u)
    us = np.array(us)
    S = np.exp(-1j * np.outer(us, g)) @ e
    Sh = np.exp(-1j * np.outer(us, g)) @ (w * e)
    proj = (S * np.conj(ph)).real; projh = (Sh * np.conj(ph)).real
    perp = (S * np.conj(ph)).imag
    se_ctrl = proj.std() / abs(p)
    se_ctrl_h = projh.std() / (abs(p) / 2)
    Gh = (w * e).sum(); bh = (Gh / (p / 2)).real
    # residual of the measured sum (along / across)
    res_al = ((G - p) * np.conj(ph)).real; res_ac = ((G - p) * np.conj(ph)).imag
    n = len(g); M = theta(T2, v['Q'], v['mu']) / np.pi
    out[key] = dict(fam=v['fam'], n=n, M=float(M), p=[p.real, p.imag], G=[G.real, G.imag], b=float(b), bh=float(bh),
                    jk={str(k): jk[k] for k in jk}, b_lo=float(b_lo), b_hi=float(b_hi), se_ctrl=float(se_ctrl),
                    se_ctrl_h=float(se_ctrl_h), ctrl_var_along_over_M=float(proj.var() / M),
                    ctrl_var_across_over_M=float(perp.var() / M), res_al=float(res_al), res_ac=float(res_ac),
                    absp=float(abs(p)))
    print('%-5s b=%.3f  jk se K=3..30: %s | ctrl se %.3f (Hann b=%.3f se %.3f) | halves %.3f %.3f | ctrl var/M along %.3f across %.3f | resid along %+.1f across %+.1f'
          % (key, b, ' '.join('%.3f' % jk[k][0] for k in jk), se_ctrl, bh, se_ctrl_h, b_lo, b_hi,
             proj.var() / M, perp.var() / M, res_al, res_ac))
json.dump(out, open(OUT + 'r2_noise.json', 'w'), indent=1)
