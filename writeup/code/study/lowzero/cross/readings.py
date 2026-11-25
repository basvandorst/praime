"""Registered readings of the cross-correlation test (round_draft.txt), one definition for the predictions, the dry run
and the test. Inputs per design: C[k][pair] (pair covariances, k over the 16 widths), the chi-values of the pairs at
2, 3, 5, 7 (psi_p = chi_1(p) chi_2(p), b_p = psi_p^2), REST[k][pair] (law, primes > 7), FULL[k][pair] (law, all primes,
finite-range form), SWAP[k][pair] (swapped term, 0 where not computed) and the per-prime tables A_T, B_T.
Per width w_k: OLS of C on Z = [1, psi_2, psi_3, psi_5, psi_7, b_2, b_3, b_5, b_7, REST, SWAP] (SWAP column only if given).
  gamma_p   = sum_k s_p(w_k) A_p(w_k) / sum_k A_p(w_k)^2               (p = 2, 3, 5, 7)
  gammaB_p  = the same for the b_p coefficients with B_p              (p = 2, 3)
  gamma_rest = mean_k coefficient of REST;  gamma_sw = mean_k coefficient of SWAP (k with sum SWAP^2 > 0)
  gamma_full = pooled slope of C on FULL over pairs and widths, centred per width
  small_w   = slope of C on FULL at w = 0.3 and 0.6 pooled (orbits up to the Heisenberg range dominate there)
  icpt      = the 16 intercepts (law: 0)
  s_p(w)    = the per-width coefficients (for the shape chi^2)
Jackknife: delete one L-function (all its pairs) at a time."""
import numpy as np

def _fit(C, Z):
    ZtZ = Z.T @ Z; Zty = Z.T @ C
    return ZtZ, Zty

def readings(Cs, PSI, REST, FULL, A, B, SWAP=None, I=None, J=None, jack=False):
    """Cs, REST, FULL, SWAP: arrays (16, npairs); PSI: (npairs, 4) values chi1(p)chi2(p) for p = 2, 3, 5, 7;
    A, B: (16, 4). Returns a dict of readings, and with jack=True also their jackknife se (needs I, J)."""
    nW, n = Cs.shape; use_sw = SWAP is not None and np.any(SWAP != 0)
    def design(k, m=None):
        cols = [np.ones(n), PSI[:, 0], PSI[:, 1], PSI[:, 2], PSI[:, 3], PSI[:, 0] ** 2, PSI[:, 1] ** 2, PSI[:, 2] ** 2, PSI[:, 3] ** 2, REST[k]]
        if use_sw: cols.append(SWAP[k])
        return np.column_stack(cols)
    mats = []
    for k in range(nW):
        Z = design(k); y = Cs[k]; x = FULL[k]
        mats.append((Z.T @ Z, Z.T @ y, np.array([n, x.sum(), y.sum(), (x * x).sum(), (x * y).sum()]), Z, y, x))
    def summarize(get):
        coef = []; S = []
        for k in range(nW):
            ZtZ, Zty, s5 = get(k); coef.append(np.linalg.lstsq(ZtZ, Zty, rcond=None)[0]); S.append(s5)
        coef = np.array(coef); S = np.array(S)
        r = {}
        for a, p in enumerate((2, 3, 5, 7)):
            r['gamma_%d' % p] = float((coef[:, 1 + a] * A[:, a]).sum() / (A[:, a] ** 2).sum())
        for a, p in enumerate((2, 3)):
            r['gammaB_%d' % p] = float((coef[:, 5 + a] * B[:, a]).sum() / (B[:, a] ** 2).sum())
        r['gamma_rest'] = float(coef[:, 9].mean())
        if use_sw:
            ok = np.array([np.sum(SWAP[k] ** 2) > 0 for k in range(nW)]); r['gamma_sw'] = float(coef[ok, 10].mean())
        nn, Sx, Sy, Sxx, Sxy = [S[:, i] for i in range(5)]
        cov = Sxy - Sx * Sy / nn; var = Sxx - Sx * Sx / nn
        r['gamma_full'] = float(cov.sum() / var.sum()); r['small_w'] = float(cov[:2].sum() / var[:2].sum())
        r['icpt'] = coef[:, 0].tolist(); r['s'] = coef[:, 1:5].tolist(); r['u'] = coef[:, 5:9].tolist()
        return r
    full = summarize(lambda k: mats[k][:3])
    if not jack: return full
    funcs = np.unique(np.concatenate([I, J])); JK = []
    masks = {}
    for f in funcs:
        m = (I == f) | (J == f)
        def get(k, m=m):
            ZtZ, Zty, s5, Z, y, x = mats[k]; Zm = Z[m]; ym = y[m]; xm = x[m]
            return ZtZ - Zm.T @ Zm, Zty - Zm.T @ ym, s5 - np.array([m.sum(), xm.sum(), ym.sum(), (xm * xm).sum(), (xm * ym).sum()])
        JK.append(summarize(get))
    nj = len(funcs); se = {}
    for key in full:
        if key in ('icpt', 's', 'u'):
            v = np.array([np.array(d[key]) for d in JK]); se[key] = np.sqrt((nj - 1) / nj * ((v - v.mean(0)) ** 2).sum(0)).tolist()
        else:
            v = np.array([d[key] for d in JK]); se[key] = float(np.sqrt((nj - 1) / nj * ((v - v.mean()) ** 2).sum()))
    # jackknife covariance of the per-width coefficients s_p(w) (for the shape chi^2)
    sv = np.array([np.array(d['s']) for d in JK]); cov_s = []
    for a in range(4):
        X = sv[:, :, a]; X = X - X.mean(0); cov_s.append(((nj - 1) / nj * X.T @ X).tolist())
    return full, se, cov_s
