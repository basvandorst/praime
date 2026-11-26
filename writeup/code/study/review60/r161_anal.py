"""Review of round 161, tasks 6-7: amplitudes on several y-ranges, jackknife over fields and over blocks of |D|,
the fixed-p two-term (Taniguchi-Thorne / BST / BTT) comparator, and the principal-form part of the spikes."""
import os, json, math
import numpy as np, mpmath as mp
H = os.path.dirname(os.path.abspath(__file__))
SC = "/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/cubic"
pred = json.load(open(SC + "/predictions_cubic.json"))   # frozen registered predictions (read only)
E = np.array(pred['edges']); h = np.diff(E); NB = 60

def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR = primes(130000)
kap = float(4 * mp.zeta(mp.mpf(1) / 3) / (5 * mp.gamma(mp.mpf(2) / 3) ** 3 * mp.zeta(mp.mpf(5) / 3))); z3 = float(mp.zeta(3))
KC = {-1: 4 * math.sqrt(3) * z3 * kap, 1: 12 * z3 * kap}   # K/C: N(X) = C X + K X^{5/6}, C = 1/(4 z3), 1/(12 z3); K = sqrt3 kap, kap

def ratio_jack(num, den, groups=None):
    """ratio estimator m = sum num / sum den per bin; delete-one jackknife over fields or over groups of fields."""
    if groups is not None:
        u = np.unique(groups); num = np.array([num[groups == g].sum(0) for g in u]); den = np.array([den[groups == g].sum(0) for g in u])
    n = len(num); N = num.sum(0); Dn = den.sum(0)
    mj = (N[None] - num) / (Dn[None] - den)
    return N / Dn, mj, n

def amp(num, den, P, sel, w, groups=None):
    m, mj, n = ratio_jack(num[:, sel], den[:, sel], groups)
    s = np.sum(P[sel] * m * w[sel]) / np.sum(P[sel] ** 2 * w[sel])
    sj = (mj * (P[sel] * w[sel])[None]).sum(1) / np.sum(P[sel] ** 2 * w[sel])
    return s, math.sqrt((n - 1) / n * np.sum((sj - sj.mean()) ** 2))

def coarse_chi2(num, den, P, sel, g=4):
    nu = num[:, sel]; de = den[:, sel]; nu = nu.reshape(len(nu), -1, g).sum(-1); de = de.reshape(len(de), -1, g).sum(-1)
    m, mj, n = ratio_jack(nu, de); dm = mj - mj.mean(0); cov = (n - 1) / n * dm.T @ dm
    Pc = P[sel].reshape(-1, g).mean(-1); r = m - Pc
    sfit = (Pc @ np.linalg.solve(cov, m)) / (Pc @ np.linalg.solve(cov, Pc)); r2 = m - sfit * Pc
    return float(r @ np.linalg.solve(cov, r)), float(r2 @ np.linalg.solve(cov, r2)), sfit

rep = {}
for sgn in (-1, 1):
    key = '%+d' % sgn; z = np.load(H + "/r161_all_%+d.npz" % sgn)
    num, den, D = z['NUM'], z['DEN'], np.abs(z['D']); P = np.array(pred['pred'][key])
    m, mj, n = ratio_jack(num, den); sig = np.sqrt((n - 1) / n * ((mj - mj.mean(0)) ** 2).sum(0)); w = 1 / sig ** 2
    R = {}
    rng = {'(0.4,3]': (E[:-1] >= 0.4 - 1e-9), '(0.2,3]': (E[:-1] >= 0.2 - 1e-9), '(0,3]': np.ones(NB, bool), '(0,0.2]': (E[1:] <= 0.2 + 1e-9),
           '(0.4,1]': (E[:-1] >= 0.4 - 1e-9) & (E[1:] <= 1 + 1e-9), '(1,2]': (E[:-1] >= 1 - 1e-9) & (E[1:] <= 2 + 1e-9), '(2,3]': (E[:-1] >= 2 - 1e-9),
           '(0.2,0.4]': (E[:-1] >= 0.2 - 1e-9) & (E[1:] <= 0.4 + 1e-9)}
    blocks = np.digitize(D, np.quantile(D, np.linspace(0, 1, 41)[1:-1]))   # 40 blocks of consecutive |D|
    for nm, sel in rng.items():
        s, se = amp(num, den, P, sel, w); _, seb = amp(num, den, P, sel, w, groups=blocks); _, sed = amp(num, den, P, sel, w, groups=D)
        R['amp ' + nm] = (round(s, 4), round(se, 4), round(seb, 4), round(sed, 4))
    # coarse chi2 on (0.4,3] (13 coarse bins of 0.2), as registered, and on (0,0.2]
    hi = rng['(0.4,3]']; lo = rng['(0,0.2]']
    R['coarse (0.4,3] chi2 pred, shape, s'] = [round(x, 2) for x in coarse_chi2(num, den, P, hi)]
    # fixed-p two-term count comparator, per field and prime
    f = pred['finiteX_factor'][key]['hold']; c_lead = (5 / 6) * KC[sgn]
    TTb = np.zeros_like(num);
    for i, A in enumerate(D):
        p = PR[(PR <= 3 * A)]; p = p[A % p != 0].astype(float); y = p / A
        x = p ** (-1 / 3); r = (5 / 6) * KC[sgn] * A ** (-1 / 6)
        T = (1 - x) / ((1 - x ** 5) * (1 + x ** 3)) * (1 + 1 / p + 1 / p ** 2)
        g = math.sqrt(A) * r * (x + x * x) * T / (1 + r * (1 + x + x * x + x ** 3) * T)
        j = np.searchsorted(E, y, side='left') - 1; ok = (j >= 0) & (j < NB)
        TTb[i] = np.bincount(j[ok], np.log(p[ok]) * g[ok], NB)
    TT = TTb.sum(0) / den.sum(0)   # same prime weights as the data
    TTauth = np.array(pred['comparators'][key]['tt_secondary_only'])
    R['c_lines (registered) vs (5/6)K/C'] = (pred['comparator_notes']['tt_secondary_only'] if False else None, round(c_lead, 6))
    R['TT exact vs authors tt comparator, max rel diff on bins'] = round(float(np.max(np.abs(TT / TTauth - 1))), 4)
    for nm, sel in (('(0,0.2]', lo), ('(0.2,0.4]', rng['(0.2,0.4]']), ('(0.4,3]', hi)):
        c2f = float(np.sum(((m - P * (f if nm == '(0,0.2]' else 1)) / sig)[sel] ** 2)); c2t = float(np.sum(((m - TT) / sig)[sel] ** 2))
        R['diag chi2 %s: lines%s, TT exact' % (nm, '*f' if nm == '(0,0.2]' else '')] = (round(c2f, 1), round(c2t, 1), int(sel.sum()))
    R['bins y<=0.4: data / lines*f / TT'] = [[round(a, 2) for a in m[:8]], [round(a, 2) for a in (P * f)[:8]], [round(a, 2) for a in TT[:8]], [round(a, 2) for a in sig[:8]]]
    # principal-form part of the spikes: n = 1, 2, 3
    for nn in range(1, 4):
        mp_ = z['NUMp'][nn - 1].sum(0) / den.sum(0)
        y0 = nn * nn / 4.0
        if sgn < 0:
            K = np.where(E[1:] > y0, (2 * np.sqrt(np.clip(E[1:] - y0, 0, None)) - 2 * np.sqrt(np.clip(E[:-1] - y0, 0, None))) / h, 0)
        else:
            K = np.where(E[:-1] < y0, (2 * np.sqrt(np.clip(y0 - E[:-1], 0, None)) - 2 * np.sqrt(np.clip(y0 - E[1:], 0, None))) / h, 0)
        sel = K > 0
        if sel.any():
            nu_ = z['NUMp'][nn - 1][:, sel]; _, mjp, nf = ratio_jack(nu_, den[:, sel])
            cf = np.sum(mp_[sel] * K[sel]) / np.sum(K[sel] ** 2); cj = (mjp * K[sel][None]).sum(1) / np.sum(K[sel] ** 2)
            R['principal-form primes n=%d: fitted coefficient of 1/sqrt|y-y0|' % nn] = (round(cf, 4), round(math.sqrt((nf - 1) / nf * np.sum((cj - cj.mean()) ** 2)), 4),
                                                                                    'shape chi2-like rms rel resid %.3f' % float(np.sqrt(np.mean(((mp_[sel] - cf * K[sel]) / (cf * K[sel])) ** 2))))
    mx = z['NUMx'].sum(0) / z['DENx'].sum(0)
    R['bins near 1/4 and 1: data, data without principal-form primes'] = [[round(a, 2) for a in m[4:8]], [round(a, 2) for a in mx[4:8]], [round(a, 2) for a in m[19:22]], [round(a, 2) for a in mx[19:22]]]
    rep[key] = R
    print(key, json.dumps(R, indent=0))
json.dump(rep, open(H + "/r161_anal.json", "w"), indent=1)
