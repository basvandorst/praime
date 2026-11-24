"""Per-prime diagonal and off-diagonal (Bogomolny-Keating 2013) changes of the number variance of Dirichlet L-zeros.
Self-contained: no imports from the repo; reads only published json summaries (r129, r120, r121, r122, r123).
usage: predict.py  (writes predict_out.txt and predict.json next to this file)"""
import os, sys, json, math, time
import numpy as np
from scipy import integrate
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bklib import (Grid, dS_diag, dsig_off, sigma2_full, sigma2_gue, C_of_r, LGRID, W, L_nodes, f_p, TWO_PI,
                   cum_sigma, R_NODES)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..'))
LZ = os.path.join(REPO, 'writeup/code/study/lowzero')
OUT = {}
PR = (2, 3, 5, 7, 11, 13)
t0 = time.time()


def idx(l):
    return int(round(l / W)) - 1


def prime_factors(n):
    n = abs(n); ps = []; p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps


G = Grid(eps_max=7.6, h=0.0025, pmax=10 ** 6)
print('grid built %.1fs' % (time.time() - t0))

# ---------- A. accuracy checks ----------
import mpmath as mp
mp.mp.dps = 30
print('\n[A] spline vs direct evaluation at off-grid eps')
from bklib import primes_upto, clog1p
Pr = primes_upto(10 ** 6).astype(float); lP = np.log(Pr)
for e in (0.0123, 1.23456, 3.3333, 6.789):
    s = mp.mpc(1, e); z = mp.zeta(s)
    Z2d = float(abs(1j * e * z) ** 2)
    w = -4 * np.sin(e * lP / 2) ** 2 * np.exp(1j * e * lP) / (Pr - 1) ** 2
    Phid = np.exp(clog1p(-w).sum())
    print('  eps %.4f  Z2 spline-direct %.2e   Phi_off spline-direct %.2e' % (e, G('Z2', e) - Z2d, abs(G('Phi', e) - Phid)))
# truncation of the prime product
Pr5 = Pr[Pr < 1e5]; lP5 = np.log(Pr5)
for e in (1.0, 5.0):
    w6 = -4 * np.sin(e * lP / 2) ** 2 * np.exp(1j * e * lP) / (Pr - 1) ** 2
    w5 = -4 * np.sin(e * lP5 / 2) ** 2 * np.exp(1j * e * lP5) / (Pr5 - 1) ** 2
    print('  eps %.1f  |Phi_off(p<1e6) - Phi_off(p<1e5)| = %.2e' % (e, abs(np.exp(clog1p(-w6).sum()) - np.exp(clog1p(-w5).sum()))))
# quadrature check of one off-diagonal integral against scipy.quad, period by period
L0 = 9.38; p0 = 2; l0 = 7.5
def integrand(r):
    e = TWO_PI * r / L0
    return 2 * (l0 - r) * (2 / L0 ** 2) * (np.exp(1j * TWO_PI * r) * G('Z2', e) * G('Phi', e) * f_p(p0, e)).real
q = sum(integrate.quad(integrand, a, min(a + 0.5, l0), limit=200, epsabs=1e-13)[0] for a in np.arange(0, l0, 0.5))
print('  dSigma2_off(p=2, L=9.38, l=7.5): Gauss pieces %.8f  scipy.quad per half-period %.8f' % (dsig_off(G, L0, 2)[idx(7.5)], q))
# GUE check against the closed form (Mehta): Sigma2 = (1/pi^2)[ln(2 pi l) + gamma + 1 - Ci(2 pi l)] + l[1 - (2/pi) Si(2 pi l)] - (1/pi^2) sin^2(pi l)...
from scipy.special import sici
def gue_closed(l):
    si, ci = sici(2 * math.pi * l)
    # exact: l - 2 int_0^l (l-r) sinc^2 ; computed by quad here as an independent check
    return l - 2 * integrate.quad(lambda x: (l - x) * (math.sin(math.pi * x) / (math.pi * x)) ** 2 if x else l, 0, l, limit=400)[0]
sg = sigma2_gue()
print('  GUE Sigma2: l=1 %.6f (quad %.6f), l=8 %.6f (quad %.6f)' % (sg[idx(1)], gue_closed(1), sg[idx(8)], gue_closed(8)))
# C(0) must vanish (R(0) = 0): check at tiny r for several bad sets
for B in ((), (2,), (3, 5), (2, 7, 13)):
    print('  C(r=1e-4) for bad set %s: %+.2e' % (B, C_of_r(G, 9.38, B, r=np.array([1e-4]))[0]))

# ---------- B. zeta check (retrospective, r129.json) ----------
print('\n[B] zeta: full BK Sigma2 (k = 1) against the measured values of round 129 (26,812 zeros, 7.6 <= ln(t/2pi) <= 8.6)')
r129 = json.load(open(os.path.join(REPO, 'writeup/data/study/r129.json')))
Ln, wn = L_nodes(1, TWO_PI * math.exp(7.6), TWO_PI * math.exp(8.6), n=12)
bk = sum(w * sigma2_full(G, L, ()) for L, w in zip(Ln, wn))
bkf = sum(w * sigma2_full(G, L, (), flip=True) for L, w in zip(Ln, wn))
bk_fixedL = sigma2_full(G, r129['L'], ())
print('  mean L of the nodes %.4f (r129: %.4f)' % (float(np.dot(Ln, wn)), r129['L']))
print('     l    data     BK   BK(fixed L)  BK(flipped sign)  Berry    GUE   data-BK  data-Berry')
zrows = []
for row in r129['rows']:
    l = row['l']; i = idx(l)
    zrows.append(dict(l=l, data=row['data'], bk=float(bk[i]), bk_fixedL=float(bk_fixedL[i]), bk_flip=float(bkf[i]), berry=row['berry'], gue=row['gue_exact']))
    print('  %4.1f  %.4f  %.4f   %.4f      %.4f          %.4f  %.4f  %+.4f  %+.4f' % (l, row['data'], bk[i], bk_fixedL[i], bkf[i], row['berry'], row['gue_exact'], row['data'] - bk[i], row['data'] - row['berry']))
d = np.array([r['data'] for r in zrows]); b = np.array([r['bk'] for r in zrows]); bb = np.array([r['berry'] for r in zrows]); bf = np.array([r['bk_flip'] for r in zrows])
print('  max |data-BK| %.4f, mean data-BK %+.4f;  max |data-Berry| %.4f, mean %+.4f;  max |data-BKflipped| %.4f' % (
    np.abs(d - b).max(), (d - b).mean(), np.abs(d - bb).max(), (d - bb).mean(), np.abs(d - bf).max()))
OUT['zeta_check'] = zrows

# ---------- C. per-prime tables ----------
def table(L, label):
    print('\n[C] per-prime change at fixed L = %.3f (%s): diag | off (exact single prime) | off (leading 1/p) | total' % (L, label))
    res = {}
    for p in PR:
        dg = np.array([dS_diag(p, L, l) for l in LGRID]); of = dsig_off(G, L, p); ofl = dsig_off(G, L, p, lead=True)
        res[p] = dict(diag=dg, off=of, off_lead=ofl)
        print('  p=%2d' % p)
        print('    l     ' + ' '.join('%7.1f' % l for l in np.arange(0.5, 8.01, 0.5)))
        for nm, arr in (('diag', dg), ('off', of), ('off1/p', ofl), ('total', dg + of)):
            print('    %-7s' % nm + ' '.join('%+7.4f' % arr[idx(l)] for l in np.arange(0.5, 8.01, 0.5)))
    return res

T938 = table(9.38, 'registered value of round 123')
OUT['table_L9.38'] = {str(p): {k: [float(v[idx(l)]) for l in np.arange(0.5, 8.01, 0.5)] for k, v in T938[p].items()} for p in PR}
for L, lab in ((9.2246, 'round 121 complex, mean L'), (10.6859, 'round 122, mean L')):
    T = {p: dict(diag=np.array([dS_diag(p, L, l) for l in LGRID]), off=dsig_off(G, L, p)) for p in PR}
    print('\n[C] L = %.4f (%s): off-diagonal part only, l = 0.5..8' % (L, lab))
    for p in PR:
        print('  p=%2d off ' % p + ' '.join('%+7.4f' % T[p]['off'][idx(l)] for l in np.arange(0.5, 8.01, 0.5)))
    OUT['table_L%.4f' % L] = {str(p): {k: [float(v[idx(l)]) for l in np.arange(0.5, 8.01, 0.5)] for k, v in T[p].items()} for p in PR}

# ---------- D. shape on the fine grid ----------
print('\n[D] shape of the off-diagonal part on the fine grid (step %.3f), L = 9.38' % W)
for p in PR:
    of = T938[p]['off']
    whole = [of[idx(l)] for l in range(1, 9)]; half = [of[idx(l + 0.5)] for l in range(0, 8)]
    # local extrema
    mx = [LGRID[i] for i in range(1, len(of) - 1) if of[i] > of[i - 1] and of[i] > of[i + 1]]
    mn = [LGRID[i] for i in range(1, len(of) - 1) if of[i] < of[i - 1] and of[i] < of[i + 1]]
    print('  p=%2d whole l 1..8: %s' % (p, ' '.join('%+.4f' % v for v in whole)))
    print('       half  l .5..7.5: %s' % ' '.join('%+.4f' % v for v in half))
    print('       maxima at l = %s' % ' '.join('%.3f' % v for v in mx[:9]))
    print('       minima at l = %s' % ' '.join('%.3f' % v for v in mn[:9]))
    sw = [of[idx(l + 0.5)] - 0.5 * (of[idx(l)] + of[idx(l + 1)]) for l in range(0, 8)]
    print('       half minus mean of neighbouring whole: %s' % ' '.join('%+.4f' % v for v in sw))
# ratio at l = 0.5 (total/diag) per prime at L=9.38
print('  l = 0.5: total/diag per prime: ' + ' '.join('%d: %.3f' % (p, (T938[p]['diag'][idx(.5)] + T938[p]['off'][idx(.5)]) / T938[p]['diag'][idx(.5)]) for p in PR))
print('  l = 1.0: total/diag per prime: ' + ' '.join('%d: %.3f' % (p, (T938[p]['diag'][idx(1)] + T938[p]['off'][idx(1)]) / T938[p]['diag'][idx(1)]) for p in PR))
OUT['fine_L9.38'] = {'l': [float(x) for x in LGRID], **{str(p): [float(x) for x in T938[p]['off']] for p in PR}}

# ---------- E. retrospective: round 123 (181 real L-functions, 100 <= gamma <= 600) ----------
print('\n[E] RETROSPECTIVE (fitted values already known): round 123 per-prime coefficients')
r123 = json.load(open(os.path.join(LZ, 'r123.json')))
r120 = json.load(open(os.path.join(LZ, 'r120.json')))
rows = r120['rows']
LS = [0.5 * i for i in range(1, 17)]
key = lambda p, l: '%d_%s' % (p, float(l))
fit = np.array([[r123['fit'][key(p, l)] for p in PR] for l in LS]); se = np.array([[r123['se'][key(p, l)] for p in PR] for l in LS])
pred_repo = np.array([[r123['pred'][key(p, l)] for p in PR] for l in LS])
# V0: our diag at each D's mean L, group mean (the repo's method) -> must equal r123 pred
pd0 = np.array([[np.mean([dS_diag(p, r['L'], l) for r in rows if r['D'] % p == 0]) for p in PR] for l in LS])
print('  check: our diagonal (group mean at each mean L) vs r123 pred: max |diff| = %.2e' % np.abs(pd0 - pred_repo).max())
# V1: group mean of the single-prime off-diagonal at each D's mean L
offcache = {}
def off_at(p, L):
    k = (p, round(L, 6))
    if k not in offcache: offcache[k] = dsig_off(G, L, p)
    return offcache[k]
po1 = np.array([[np.mean([off_at(p, r['L'])[idx(l)] for r in rows if r['D'] % p == 0]) for p in PR] for l in LS])
# V2: regression-matched: predicted Sigma2 per D (all its primes, exact product, L-averaged over 100..600), same OLS as r123
X = np.array([[1] + [int(r['D'] % p == 0) for p in PR] for r in rows], float); XtXi = np.linalg.inv(X.T @ X)
def predicted_sigma(Ds, t1, t2, diag_only=False):
    Y = []
    for D in Ds:
        q = abs(D); B = tuple(prime_factors(D)); Ln, wn = L_nodes(q, t1, t2, n=8)
        s = 0
        for L, w in zip(Ln, wn):
            if diag_only:
                s = s + w * (sigma2_full(G, L, ()) + sum(np.array([dS_diag(p, L, l) for l in LGRID]) for p in B))
            else:
                s = s + w * sigma2_full(G, L, B)
        Y.append(s)
    return np.array(Y)
Ds = [r['D'] for r in rows]
Yfull = predicted_sigma(Ds, 100, 600); Ydiag = predicted_sigma(Ds, 100, 600, diag_only=True)
def regress(Y):
    out = np.zeros((len(LS), len(PR)))
    for j, l in enumerate(LS):
        c = XtXi @ X.T @ Y[:, idx(l)]; out[j] = c[1:]
    return out
pv2_full = regress(Yfull); pv2_diag = regress(Ydiag)
def chi2(pred):
    z = (fit - pred) / se
    return (z ** 2).sum(), (z ** 2).sum(axis=0), (z ** 2).sum(axis=1)
variants = [('diag only (r123 pred)', pred_repo), ('diag + off, group mean at mean L', pd0 + po1),
            ('diag only, regression-matched, L-averaged', pv2_diag), ('diag + off, regression-matched, L-averaged (full BK)', pv2_full)]
OUT['retro'] = {}
for nm, pv in variants:
    tot, perp, perl = chi2(pv)
    sl = float((fit * pv).sum() / (pv * pv).sum())
    print('  %-55s chi2 = %7.1f on 96 | per prime %s | slope %.3f' % (nm, tot, ' '.join('%d:%.1f' % (p, c) for p, c in zip(PR, perp)), sl))
    print('  %-55s chi2 per l (6 dof): %s' % ('', ' '.join('%.0f' % c for c in perl)))
    OUT['retro'][nm] = dict(chi2=float(tot), per_prime={str(p): float(c) for p, c in zip(PR, perp)}, per_l=[float(c) for c in perl], slope=sl, pred=pv.tolist())
# l = 0.5 residuals, and the slope at l = 0.5
j = 0
print('  l = 0.5: fit         ' + ' '.join('%+.4f' % v for v in fit[j]) + '   (se ' + ' '.join('%.4f' % v for v in se[j]) + ')')
for nm, pv in variants:
    print('  l = 0.5: %-14s' % nm[:14] + ' '.join('%+.4f' % v for v in pv[j]) + '  resid ' + ' '.join('%+.4f' % v for v in (fit[j] - pv[j])) + '  slope fit/pred %.2f' % ((fit[j] @ pv[j]) / (pv[j] @ pv[j])))
# residual of the diagonal law against the predicted off-diagonal part, slope through the origin
res_d = fit - pred_repo; offp = pv2_full - pv2_diag
g = float((res_d * offp).sum() / (offp * offp).sum()); gw = float(((res_d * offp) / se ** 2).sum() / ((offp ** 2) / se ** 2).sum())
print('  residual(fit - diag) against predicted off-diagonal part over 96 cells: slope %.3f (se-weighted %.3f), corr %.3f' % (g, gw, np.corrcoef(res_d.ravel(), offp.ravel())[0, 1]))
OUT['retro']['resid_vs_off'] = dict(slope=g, slope_w=gw, corr=float(np.corrcoef(res_d.ravel(), offp.ravel())[0, 1]))
# the review's example: p = 2 at l = 7 and 7.5
for l in (7.0, 7.5):
    j = LS.index(l)
    print('  p=2, l=%.1f: fit - diag %+.4f (se %.4f); predicted off (regression-matched) %+.4f; group-mean off %+.4f' % (l, res_d[j, 0], se[j, 0], offp[j, 0], po1[j, 0]))
print('  per prime, residual (fit - diag) and predicted off, l = 0.5..8:')
for i, p in enumerate(PR):
    print('   p=%2d resid ' % p + ' '.join('%+.4f' % v for v in res_d[:, i]))
    print('        off   ' + ' '.join('%+.4f' % v for v in offp[:, i]))
    print('        se    ' + ' '.join('%.4f' % v for v in se[:, i]))
# absolute check at l = 1 (published Sigma2(1) per D)
s1 = np.array([r['sigma2_1'] for r in rows])
print('  absolute Sigma2(1): mean measured %.4f, mean BK %.4f, mean diff %+.4f, corr %.3f, rms diff %.4f' % (
    s1.mean(), Yfull[:, idx(1)].mean(), (s1 - Yfull[:, idx(1)]).mean(), np.corrcoef(s1, Yfull[:, idx(1)])[0, 1], np.sqrt(((s1 - Yfull[:, idx(1)]) ** 2).mean())))
OUT['abs_r120_l1'] = dict(meas=float(s1.mean()), bk=float(Yfull[:, idx(1)].mean()))
# rounds 121 (complex, both sides 100..500) and 122 (real, 900..1500) at l = 1, absolute and beta
for name, fn, t1, t2, keyN in (('r121 complex', 'r121.json', 100, 500, 'N'), ('r122 real 900-1500', 'r122.json', 900, 1500, 'D')):
    rr = json.load(open(os.path.join(LZ, fn)))['rows']
    Dn = [r[keyN] for r in rr]
    Yf = predicted_sigma(Dn, t1, t2); Yd = predicted_sigma(Dn, t1, t2, diag_only=True)
    m = np.array([r['sigma2_1'] for r in rr])
    xs = np.array([r['x'] for r in rr])
    A = np.column_stack([np.ones(len(m)), xs])
    bm = np.linalg.lstsq(A, m, rcond=None)[0][1]; bf = np.linalg.lstsq(A, Yf[:, idx(1)], rcond=None)[0][1]; bd = np.linalg.lstsq(A, Yd[:, idx(1)], rcond=None)[0][1]
    print('  %s l=1: mean measured %.4f, BK %.4f (diff %+.4f); beta measured %.3f, beta of BK-predicted values %.3f, of diag-predicted %.3f' % (
        name, m.mean(), Yf[:, idx(1)].mean(), (m - Yf[:, idx(1)]).mean(), bm, bf, bd))
    OUT['abs_' + name] = dict(meas=float(m.mean()), bk=float(Yf[:, idx(1)].mean()), beta_meas=float(bm), beta_bk=float(bf), beta_diag=float(bd))

# ---------- F. predictions for the fresh tests ----------
print('\n[F] PREDICTIONS for the fresh tests (regression-matched, L-averaged, exact product over bad primes)')
def fund_discs(lo, hi):
    out = []
    for D in list(range(-hi, -lo + 1)) + list(range(lo, hi + 1)):
        if D in (0, 1): continue
        if D % 4 == 1 and all(D % (p * p) for p in range(2, int(abs(D) ** .5) + 1)): out.append(D)
        elif D % 4 == 0:
            m = D // 4
            if m % 4 in (2, 3) and all(m % (p * p) for p in range(2, int(abs(m) ** .5) + 2)): out.append(D)
    return sorted(out, key=lambda x: (abs(x), x))
# sanity: count in 100 <= |D| <= 400 must be 181
print('  fundamental discriminants with 100 <= |D| <= 400: %d (round 120: 181)' % len(fund_discs(100, 400)))
def design_and_pred(Ds_, primes_, t1, t2):
    Xd = np.array([[1] + [int(D % p == 0) for p in primes_] for D in Ds_], float); XtXi_ = np.linalg.inv(Xd.T @ Xd)
    Yf = predicted_sigma(Ds_, t1, t2); Yd = predicted_sigma(Ds_, t1, t2, diag_only=True)
    cf = np.array([(XtXi_ @ Xd.T @ Yf[:, i])[1:] for i in range(len(LGRID))]); cd = np.array([(XtXi_ @ Xd.T @ Yd[:, i])[1:] for i in range(len(LGRID))])
    return cf, cd, Yf, Yd, Xd
FRESH = {}
# Test A: 88 complex characters of round 121 (conductors from r121.json), indicators of 2, 3, 5, 7
rr = json.load(open(os.path.join(LZ, 'r121.json')))['rows']
NsA = [r['N'] for r in rr]
cfA, cdA, YfA, YdA, XA = design_and_pred(NsA, (2, 3, 5, 7), 100, 500)
print('  Test A (88 complex, 100 <= |gamma| <= 500): group sizes %s' % [int(XA[:, i + 1].sum()) for i in range(4)])
FRESH['A'] = dict(primes=[2, 3, 5, 7], n=len(NsA), full=cfA.tolist(), diag=cdA.tolist())
# Test B: new real characters 400 < |D| <= 700, zeros 100 <= gamma <= 600, indicators of 2..13
DsB = fund_discs(401, 700)
cfB, cdB, YfB, YdB, XB = design_and_pred(DsB, PR, 100, 600)
LB = [float(np.dot(*L_nodes(abs(D), 100, 600, 8))) for D in DsB]
print('  Test B (real, 400 < |D| <= 700, 100 <= gamma <= 600): n = %d, mean L %.3f, group sizes %s' % (len(DsB), np.mean(LB), [int(XB[:, i + 1].sum()) for i in range(6)]))
FRESH['B'] = dict(primes=list(PR), n=len(DsB), D=DsB, meanL=float(np.mean(LB)), full=cfB.tolist(), diag=cdB.tolist())
# Test C (fallback): the 181 discriminants of round 120 at the new band 600 <= gamma <= 1000
cfC, cdC, YfC, YdC, XC = design_and_pred(Ds, PR, 600, 1000)
LC = [float(np.dot(*L_nodes(abs(D), 600, 1000, 8))) for D in Ds]
print('  Test C (round-120 discriminants, 600 <= gamma <= 1000): n = %d, mean L %.3f' % (len(Ds), np.mean(LC)))
FRESH['C'] = dict(primes=list(PR), n=len(Ds), meanL=float(np.mean(LC)), full=cfC.tolist(), diag=cdC.tolist())
for nm, cf, cd, prs in (('A', cfA, cdA, (2, 3, 5, 7)), ('B', cfB, cdB, PR), ('C', cfC, cdC, PR)):
    print('\n  Test %s: predicted coefficient per prime, full BK (diag + off) / diag only, l = 0.5..8' % nm)
    print('    l    ' + ' '.join('%15s' % ('p=%d' % p) for p in prs))
    for l in np.arange(0.5, 8.01, 0.5):
        i = idx(l)
        print('   %4.1f  ' % l + ' '.join('%+.4f/%+.4f' % (cf[i, k], cd[i, k]) for k in range(len(prs))))
    # off-diagonal part summary
    offp_ = cf - cd
    half = [idx(l) for l in np.arange(4.5, 7.51, 1.0)]; whole = [idx(l) for l in np.arange(4.0, 8.01, 1.0)]
    print('    pooled over primes (sum of coefficients) of the off part: mean over half l 4.5..7.5 = %+.4f, mean over whole l 4..8 = %+.4f' % (
        offp_[half].sum(1).mean(), offp_[whole].sum(1).mean()))
    print('    l = 0.5 slope (full vs diag) %.3f' % ((cf[idx(.5)] @ cd[idx(.5)]) / (cd[idx(.5)] @ cd[idx(.5)])))
    # registered statistics, their values under the two laws
    LSr = np.arange(0.5, 8.01, 0.5)
    o = np.array([cf[idx(l)] - cd[idx(l)] for l in LSr])
    print('    p=2 off part (full - diag) at l = 6, 6.5, 7, 7.5, 8: ' + ' '.join('%+.4f' % o[list(LSr).index(l), 0] for l in (6, 6.5, 7, 7.5, 8)))
    FRESH[nm]['l05_ratio'] = float((cf[idx(.5)] @ cd[idx(.5)]) / (cd[idx(.5)] @ cd[idx(.5)]))
    FRESH[nm]['p2_off_6_to_8'] = [float(o[list(LSr).index(l), 0]) for l in (6, 6.5, 7, 7.5, 8)]
# r123 design: predicted l = 0.5 ratio full/diag (regression-matched)
print('\n  round 123 design: predicted l = 0.5 ratio full/diag %.3f; fitted/diag %.3f' % (
    (pv2_full[0] @ pv2_diag[0]) / (pv2_diag[0] @ pv2_diag[0]), (fit[0] @ pv2_diag[0]) / (pv2_diag[0] @ pv2_diag[0])))
OUT['r123_l05_ratio_pred'] = float((pv2_full[0] @ pv2_diag[0]) / (pv2_diag[0] @ pv2_diag[0]))
# per-L-function predictions for the runner, l = 0.125 .. 8 step 0.125: full BK and diag-only law (base + Berry drops)
sel = [idx(l) for l in np.arange(0.125, 8.001, 0.125)]
PERL = {'l': [float(LGRID[i]) for i in sel],
        'A': {'height': [100, 500], 'sides': 'both', 'N': NsA, 'full': YfA[:, sel].tolist(), 'diag': YdA[:, sel].tolist()},
        'B': {'height': [100, 600], 'D': DsB, 'full': YfB[:, sel].tolist(), 'diag': YdB[:, sel].tolist()},
        'C': {'height': [600, 1000], 'D': Ds, 'full': YfC[:, sel].tolist(), 'diag': YdC[:, sel].tolist()}}
json.dump(PERL, open(os.path.join(HERE, 'predictions_per_lfunction.json'), 'w'))
json.dump({'out': OUT, 'fresh': FRESH, 'lgrid': LGRID.tolist()}, open(os.path.join(HERE, 'predict.json'), 'w'))
print('\ndone %.1fs' % (time.time() - t0))
