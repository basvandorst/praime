"""Round 104 follow-up (not pre-registered): how much of the first zero is fixed by the exact identity
S = sum_rho 1/|rho|^2 = 2 L'/L(1, chi) + ln(q/pi) + psi((1+a)/2)? Compare gamma_1 with the prediction obtained when the
first pair takes the excess of S over the median S of discriminants of similar size: 2/(1/4 + g^2) = S - S_med + 2/(1/4 + gmed^2).
usage: r104b.py"""
import json, math, os
import numpy as np, cypari2
from scipy.stats import spearmanr
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r104.json')))
rows = d['rows']
for r in rows:
    L = pari.lfuncreate(r['D'])
    r['dlog'] = float(pari.lfun(L, 1, 1)) / r['L1']
    a = 1 if r['D'] < 0 else 0
    r['S'] = 2 * r['dlog'] + math.log(abs(r['D']) / math.pi) + float(pari.psi((1 + a) / 2.0))
out = {}
for name, sel in (('negative', lambda r: r['D'] < 0), ('positive', lambda r: r['D'] > 0)):
    R = [r for r in rows if sel(r)]
    rho = spearmanr([r['S'] for r in R], [r['gamma1'] for r in R])[0]
    # prediction: compare within bins of ln|D|
    lnq = np.array([math.log(abs(r['D'])) for r in R]); S = np.array([r['S'] for r in R]); g = np.array([r['gamma1'] for r in R])
    pred = np.empty(len(R))
    for i in range(len(R)):
        m = np.abs(lnq - lnq[i]) < 0.25
        Smed = np.median(S[m]); gmed = np.median(g[m])
        lhs = S[i] - Smed + 2 / (0.25 + gmed ** 2)
        pred[i] = math.sqrt(max(2 / lhs - 0.25, 0)) if lhs > 0 else np.nan
    ok = ~np.isnan(pred)
    rho2 = spearmanr(pred[ok], g[ok])[0]
    lo = np.argsort(g)[:200]
    out[name] = {'spearman_S_gamma1': rho, 'spearman_pred_gamma1': rho2,
                 'median_abs_err_lowest200': float(np.nanmedian(np.abs(pred[lo] - g[lo]))), 'median_gamma1_lowest200': float(np.median(g[lo]))}
    print(name, out[name])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r104b.json'), 'w'), indent=1)
