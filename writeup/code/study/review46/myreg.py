"""Own regression readings (definitions as registered in round141.txt) and delete-one-function jackknife."""
import numpy as np


def fit_all(Cp, PSI, REST, FULL, SW, A, B, extra=None):
    """Cp, REST, FULL, SW: (16, npairs). PSI (npairs, 4). Returns dict of readings and per-width coefs."""
    nW, n = Cp.shape
    use_sw = SW is not None and np.any(SW != 0)
    coefs = []
    for k in range(nW):
        cols = [np.ones(n)] + [PSI[:, a] for a in range(4)] + [PSI[:, a] ** 2 for a in range(4)] + [REST[k]]
        if use_sw: cols.append(SW[k])
        if extra is not None: cols += [e[k] for e in extra]
        Z = np.column_stack(cols)
        coefs.append(np.linalg.lstsq(Z, Cp[k], rcond=None)[0])
    coefs = np.array(coefs)
    r = {}
    for a, p in enumerate((2, 3, 5, 7)):
        r['gamma_%d' % p] = (coefs[:, 1 + a] * A[:, a]).sum() / (A[:, a] ** 2).sum()
    for a, p in enumerate((2, 3)):
        r['gammaB_%d' % p] = (coefs[:, 5 + a] * B[:, a]).sum() / (B[:, a] ** 2).sum()
    r['gamma_rest'] = coefs[:, 9].mean()
    if use_sw: r['gamma_sw'] = coefs[:, 10].mean()
    xc = FULL - FULL.mean(1, keepdims=True); yc = Cp - Cp.mean(1, keepdims=True)
    cv = (xc * yc).sum(1); vv = (xc * xc).sum(1)
    r['gamma_full'] = cv.sum() / vv.sum(); r['small_w'] = cv[:2].sum() / vv[:2].sum()
    r['icpt_mean'] = coefs[:, 0].mean()
    return r, coefs


def jackknife(Cp, PSI, REST, FULL, SW, A, B, I, J, extra=None):
    full, coefs = fit_all(Cp, PSI, REST, FULL, SW, A, B, extra)
    funcs = np.unique(np.concatenate([I, J])); reps = []; creps = []
    for f in funcs:
        m = ~((I == f) | (J == f))
        r, c = fit_all(Cp[:, m], PSI[m], REST[:, m], FULL[:, m], None if SW is None else SW[:, m], A, B,
                       None if extra is None else [e[:, m] for e in extra])
        reps.append(r); creps.append(c)
    nj = len(funcs); se = {}
    for key in full:
        v = np.array([d[key] for d in reps]); se[key] = np.sqrt((nj - 1) / nj * ((v - v.mean()) ** 2).sum())
    cr = np.array(creps); cse = np.sqrt((nj - 1) / nj * ((cr - cr.mean(0)) ** 2).sum(0))
    return full, se, coefs, cse, cr
