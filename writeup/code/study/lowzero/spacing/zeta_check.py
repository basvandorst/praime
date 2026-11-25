"""Validation on zeta (task 2): Var(s) and <r~> of zeta's zeros against the level predictions (level.py; second-order
deformation formula with the exact two-point function), BBLM's 1/L^2 (+Q) formula, and the sine kernel.
Data: the first 40,000 zeros (writeup/data/riemann/zeros_40000_flint.txt), and Odlyzko's 10^4 zeros near 10^12, 10^21, 10^22.
Unfolding u = (t/2pi) ln(t/(2 pi e)) + 7/8 (+1/(48 pi t)); spacings s_n = u_{n+1} - u_n; r~_n = min(s_n, s_{n+1})/max(...).
Errors: blocks of 200 consecutive spacings (block bootstrap, 2000 resamples)."""
import sys, os, math, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import level
R = os.path.join(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..')), 'data', 'riemann') + '/'
rng = np.random.default_rng(7)

def unfold(t):
    return t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 7 / 8 + 1 / (48 * math.pi * t)

def stats(s):
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return np.var(s), np.mean(r)

def boot(s, nb=200, B=2000):
    nblk = len(s) // nb; blocks = [s[i * nb:(i + 1) * nb] for i in range(nblk)]
    est = []
    for _ in range(B):
        idx = rng.integers(0, nblk, nblk); x = np.concatenate([blocks[i] for i in idx]); est.append(stats(x))
    return np.array(est).std(0)

def pred_avg(Ls):
    """predictions averaged over the spacings' L values (the level is smooth in L: average over 7 quantile nodes)"""
    qs = np.quantile(Ls, [1 / 14 + k / 7 for k in range(7)])
    v = np.array([level.level(L) for L in qs]).mean(0)
    lam, Q, V0 = 1.573151, 2.315846, level.V0
    bb = np.mean(V0 - (2 * lam / Ls ** 2) * (1 + V0)); bq = np.mean(V0 - (2 * lam / Ls ** 2) * (1 + V0) + 6 * Q * (1 + V0) / Ls ** 3)
    return level.V0 + v[0], level.R0 + v[1], bb, bq

out = []
t = np.loadtxt(R + 'zeros_40000_flint.txt'); u = unfold(t); s = np.diff(u); mid = 0.5 * (t[1:] + t[:-1]); Lm = np.log(mid / (2 * math.pi))
print('zeta, first 40,000 zeros: bins in L = ln(t/2pi)')
print('  L range        n      Var(s)  se       pred    BBLM    BBLM+Q |  <r~>    se      pred    (sine 0.17999 / 0.59975)')
for lo, hi in ((4.0, 5.5), (5.5, 6.5), (6.5, 7.5), (7.5, 8.0), (8.0, 8.6), (7.6, 8.6)):
    sel = (Lm >= lo) & (Lm < hi); x = s[sel]
    v, r = stats(x); se = boot(x)
    pv, pr, bb, bq = pred_avg(Lm[sel])
    out.append(dict(set='zeta40k', L=[lo, hi], n=int(sel.sum()), meanL=float(Lm[sel].mean()), var=v, var_se=se[0], r=r, r_se=se[1], pred_var=pv, pred_r=pr, bblm=bb, bblmq=bq))
    print('  %4.1f-%4.1f  %6d  %.5f %.5f  %.5f %.5f %.5f | %.5f %.5f  %.5f' % (lo, hi, len(x), v, se[0], pv, bb, bq, r, se[1], pr))
print('\nOdlyzko sets (10^4 zeros each)')
for f, off in (('odlyzko_zeros3.txt', 267653395647), ('odlyzko_zeros4.txt', 144176897509546973000), ('odlyzko_zeros5.txt', 1370919909931995300000)):
    vals = [ln.strip() for ln in open(R + f) if re.fullmatch(r'[0-9]+\.[0-9]+', ln.strip())]
    d = np.array([float(x) for x in vals])                 # gamma - off
    g0 = float(off); L = math.log(g0 / (2 * math.pi))
    s = np.diff(d) * L / (2 * math.pi)                      # local unfolding: density L/2pi (constant across the window)
    v, r = stats(s); se = boot(s, nb=100)
    dv, dr = level.level(L); lam, Q, V0 = 1.573151, 2.315846, level.V0
    out.append(dict(set=f, L=L, n=len(s), var=v, var_se=se[0], r=r, r_se=se[1], pred_var=V0 + dv, pred_r=level.R0 + dr))
    print('  %s L=%.2f n=%d  Var %.5f +- %.5f (pred %.5f, BBLM %.5f)   <r~> %.5f +- %.5f (pred %.5f)  mean s %.5f'
          % (f, L, len(s), v, se[0], V0 + dv, V0 - (2 * lam / L ** 2) * (1 + V0), r, se[1], level.R0 + dr, s.mean()))
json.dump(out, open('zeta_check.json', 'w'), indent=1)
