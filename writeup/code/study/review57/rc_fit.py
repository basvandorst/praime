"""Own refit (review of round 155): time-domain weighted least squares on per-cell increments,
x in (X[1], X[4716]] = (2.0002e9, 1.9995e10], using only my own recount rc_full.npz (classes mod 420).
Model:  e_k / sqrt(X_k) = c0 * r_k / sqrt(X_k) + sum_chi beta_chi * T_chi,k / sqrt(X_k) + nuisance + noise,
  e_k = twin weight in cell k - S(h) dx_k;  r_k = theta-weight in cell k - dx_k;  T_chi,k = sum_p chi(p) ln p over the cell.
Noise: white in cells (random-walk error in u), variance prop. to (ln x)^2.  Errors: OLS and block jackknife.
usage: rc_fit.py rc_full.npz
"""
import sys, math, itertools
import numpy as np

D = np.load(sys.argv[1])
X = D['X'].astype(np.float64)
K = np.arange(int(D['i0']) + 1, int(D['i1']) + 1)   # cells (X[k-1], X[k]]
Xk = X[K]; dx = X[K] - X[K - 1]; u = np.log(Xk)
sq = np.sqrt(Xk); w = 1.0 / u ** 2
M = 420
thc = D['thc'][K]; twc = D['twc'][K]; cuc = D['cuc'][K]
th = D['th'][K]
assert np.allclose(thc.sum(1), th, rtol=1e-12)
C2 = 0.6601618158468695739
S2 = 2 * C2
units = [c for c in range(M) if math.gcd(c, M) == 1]

# ---------- characters mod 420 = 4*3*5*7 ----------
def dlog_table(p, g):
    t = {}; x = 1
    for k in range(p - 1):
        t[x] = k; x = x * g % p
    return t
comp = [(4, None), (3, 2), (5, 2), (7, 3)]
tabs = {p: dlog_table(p, g) for p, g in comp if g}
def comp_val(p, j, c):
    if p == 4:
        if c % 2 == 0: return 0
        return 1 if j == 0 else (1 if c % 4 == 1 else -1)
    if c % p == 0: return 0
    return np.exp(2j * np.pi * j * tabs[p][c % p] / (p - 1))
chars = []
for j4, j3, j5, j7 in itertools.product(range(2), range(2), range(4), range(6)):
    vals = np.array([comp_val(4, j4, c) * comp_val(3, j3, c) * comp_val(5, j5, c) * comp_val(7, j7, c) for c in range(M)], dtype=complex)
    cond = (4 if j4 else 1) * (3 if j3 else 1) * (5 if j5 else 1) * (7 if j7 else 1)
    par = int(round((vals[M - 1]).real))
    chars.append(dict(j=(j4, j3, j5, j7), v=vals, q=cond, even=(par == 1)))

def mu(q):
    r = 1
    for p in [2, 3, 5, 7]:
        if q % (p * p) == 0: return 0
        if q % p == 0: r = -r
    return r
def S_h(h):
    s = S2
    for p in [3, 5, 7, 11, 13]:
        if h % p == 0: s *= (p - 1) / (p - 2)
    return s

def beta_enum(ch, h):
    """S(h) * sum_i mean over admissible n mod 420 of conj chi(n + b_i), b = (0, h)."""
    adm = [n for n in range(M) if math.gcd(n * (n + h), M) == 1]
    c1 = np.mean([np.conj(ch['v'][n % M]) for n in adm]); c2 = np.mean([np.conj(ch['v'][(n + h) % M]) for n in adm])
    return S_h(h) * (c1 + c2)
def beta_closed(ch, h):
    q = ch['q']
    if not ch['even'] or q % 2 == 0 or mu(q) == 0 or math.gcd(q, h) != 1: return 0j
    prod = 1
    for p in [3, 5, 7]:
        if q % p == 0: prod *= (p - 2)
    j4, j3, j5, j7 = ch['j']
    val = 1
    for (p, jj) in ((4, j4), (3, j3), (5, j5), (7, j7)):
        if jj: val *= comp_val(p, jj, h)
    return 2 * S_h(h) * mu(q) * np.conj(val) / prod
for h in (2, 4):
    md = max(abs(beta_enum(ch, h) - beta_closed(ch, h)) for ch in chars if ch['q'] > 1)
    print('h=%d: closed form vs my enumeration of admissible classes, max |diff| over 95 characters = %.2e' % (h, md))

def name(ch):
    return 'q=%d j=%s %s' % (ch['q'], ch['j'], 'even' if ch['even'] else 'odd')

# ---------- regressors ----------
V = np.array([[ch['v'][c] for c in units] for ch in chars])         # 96 x 48
TH = thc[:, units] @ V.T                                             # cells x 96 complex: theta-increment of chi
r = (th - dx) / sq
nonpr = [i for i, ch in enumerate(chars) if ch['q'] > 1]
# real basis: for real chi one column; for a complex pair (chi, conj chi) keep the one with smaller index: columns 2Re, -2Im
cols, colinfo = [], []
seen = set()
for i in nonpr:
    ch = chars[i]
    if i in seen: continue
    vconj = np.conj(ch['v'])
    jj = [k for k in nonpr if np.allclose(chars[k]['v'], vconj)][0]
    seen.add(i); seen.add(jj)
    if jj == i:
        cols.append(TH[:, i].real / sq); colinfo.append((i, 're', True))
    else:
        cols.append(2 * TH[:, i].real / sq); colinfo.append((i, 're', False))
        cols.append(-2 * TH[:, i].imag / sq); colinfo.append((i, 'im', False))
CH = np.array(cols).T
uc = (u - u.mean()) / (u.max() - u.min())
NUI = np.column_stack([np.ones_like(u), uc, uc ** 2, np.exp(-(u - u[0]) / 2)])
print('cells', len(K), 'x from %.4e to %.4e' % (X[K[0] - 1], X[K[-1]]), '; character regressors', CH.shape[1])

def wls(y, Z, w, blocks=40):
    sw = np.sqrt(w)
    A = Z * sw[:, None]; b = y * sw
    G = np.linalg.inv(A.T @ A)
    beta = G @ (A.T @ b)
    res = b - A @ beta
    s2 = res @ res / (len(b) - Z.shape[1])
    cov = s2 * G
    # block jackknife (contiguous blocks in u)
    n = len(y); edges = np.linspace(0, n, blocks + 1).astype(int)
    bs = []
    for g in range(blocks):
        m = np.ones(n, bool); m[edges[g]:edges[g + 1]] = False
        Ag, bg = A[m], b[m]
        bs.append(np.linalg.solve(Ag.T @ Ag, Ag.T @ bg))
    bs = np.array(bs)
    covjk = (blocks - 1) / blocks * ((bs - bs.mean(0)).T @ (bs - bs.mean(0)))
    return beta, cov, covjk, res / sw, s2

def pred_vector(h):
    out = []
    for (i, part, isreal) in colinfo:
        b = beta_enum(chars[i], h)
        out.append(b.real if part == 're' else b.imag)
    return np.array(out)

def report(label, y, h, w=w, extra=''):
    Z = np.column_stack([r, CH, NUI])
    beta, cov, covjk, res, s2 = wls(y, Z, w)
    se, sej = np.sqrt(np.diag(cov)), np.sqrt(np.diag(covjk))
    pv = pred_vector(h)
    print('\n=== %s (h=%d) full model: zeta + all 95 characters mod 420 free + 4 nuisance; WLS on %d cells' % (label, h, len(y)))
    print('  zeta  %.4f +- %.4f (OLS)  +- %.4f (jackknife 40 blocks)   predicted 2S(h) = %.4f' % (beta[0], se[0], sej[0], 2 * S_h(h)))
    want = {(0, 0, 2, 0): 'chi_5 (5.1 real even)', (0, 1, 0, 0): 'chi_-3', (1, 0, 0, 0): 'chi_-4', (0, 0, 0, 2): 'q7 cubic', (0, 1, 1, 0): 'q15 quartic', (0, 1, 0, 3): 'chi_21 real', (0, 0, 1, 0): 'q5 odd pair', (0, 0, 0, 3): 'chi_-7'}
    for kk, (i, part, isreal) in enumerate(colinfo):
        j = chars[i]['j']
        if j in want:
            print('  %-22s %s  %+.4f +- %.4f (OLS) +- %.4f (jk)   predicted %+.4f' % (want[j], part, beta[1 + kk], se[1 + kk], sej[1 + kk], pv[kk]))
    # joint chi2 of all 95 character coefficients against the prediction, and of the null set alone
    b = beta[1:1 + CH.shape[1]]; Cb = cov[1:1 + CH.shape[1], 1:1 + CH.shape[1]]; Cj = covjk[1:1 + CH.shape[1], 1:1 + CH.shape[1]]
    d = b - pv
    chi2_all = d @ np.linalg.solve(Cb, d)
    null = np.array([abs(pv[k]) < 1e-12 for k in range(len(pv))])
    dn = b[null]
    chi2_null = dn @ np.linalg.solve(Cb[np.ix_(null, null)], dn)
    _, _, covjk200, _, _ = wls(y, Z, w, blocks=235)
    Cj2 = covjk200[1:1 + CH.shape[1], 1:1 + CH.shape[1]]
    chi2_null_j = dn @ np.linalg.solve(Cj2[np.ix_(null, null)], dn)
    zz = [(d[k] / math.sqrt(Cb[k, k]), name(chars[colinfo[k][0]]), colinfo[k][1]) for k in range(len(d)) if null[k]]
    zz.sort(key=lambda t: -abs(t[0]))
    print('  largest null z:', '; '.join('%s %s %+.2f' % (n_, p_, z_) for z_, n_, p_ in zz[:8]))
    import scipy.stats as st
    print('  null chi2 p-value (OLS): %.3f' % st.chi2.sf(chi2_null, null.sum()))
    pres = ~null
    print('  all 95 character coefs vs prediction: chi2 %.1f on %d dof (OLS cov)' % (chi2_all, len(d)))
    print('  null set (%d coefs predicted 0): chi2 %.1f (OLS cov), %.1f (jackknife cov, 235 blocks of ~20 cells)' % (null.sum(), chi2_null, chi2_null_j))
    # loudness: one scale lambda on the predicted combination (present chars, conductors 5,7,15,21,35,105)
    for sel_name, condset in [('q in {5,7,15,21}', {5, 7, 15, 21}), ('q in {35,105} (not in the registered fit)', {35, 105}), ('all present q | 105', {5, 7, 15, 21, 35, 105})]:
        sel = np.array([pres[k] and chars[colinfo[k][0]]['q'] in condset for k in range(len(pv))])
        comb = CH[:, sel] @ pv[sel]
        rest = [k for k in range(CH.shape[1]) if not sel[k]]
        Z2 = np.column_stack([r, comb, CH[:, rest], NUI])
        b2, c2, cj2, _, _ = wls(y, Z2, w)
        print('  loudness %-42s %.4f +- %.4f (OLS) +- %.4f (jk)  [zeta %.4f]' % (sel_name, b2[1], math.sqrt(c2[1, 1]), math.sqrt(cj2[1, 1]), b2[0]))
    # residual whiteness
    z = res * np.sqrt(w)
    z = z / z.std()
    ac = [np.corrcoef(z[:-L], z[L:])[0, 1] for L in (1, 2, 3, 5, 10, 50)]
    print('  residual autocorrelation lags 1,2,3,5,10,50:', ' '.join('%+.3f' % a for a in ac), ' (1/sqrt(n) = %.3f)' % (1 / math.sqrt(len(z))))
    # variance ratio of residual block sums (random walk check): sums over blocks of m cells
    for m in (8, 64, 512):
        nb = len(z) // m
        bsums = z[:nb * m].reshape(nb, m).sum(1)
        print('    block-sum variance ratio m=%d: %.3f' % (m, bsums.var() / m), end='')
    print()
    return beta, cov, covjk, res

e2 = (twc[:, units].sum(1) - S2 * dx) / sq
e4 = (cuc[:, units].sum(1) - S_h(4) * dx) / sq
tw_total = D['P2'][K]
assert np.allclose(twc.sum(1), tw_total, rtol=1e-12)
b2, cov2, covj2, res2 = report('twins', e2, 2)
b4, cov4, covj4, res4 = report('cousins', e4, 4)
np.savez(sys.argv[2] if len(sys.argv) > 2 else 'rc_fit_out.npz', res2=res2, res4=res4, r=r, CH=CH, NUI=NUI, w=w, u=u, b2=b2, b4=b4)
