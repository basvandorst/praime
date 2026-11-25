"""tzlib.py -- series on the ln-x grid and the band-limited regression of a pattern error on the prime errors.

Series (all divided by sqrt(x), as functions of u = ln x on the grid edges X_j):
  Y_pat(u)   = (sum_{n<=x} prod Lambda-weights - S x) / sqrt(x)        (primes only: weights ln p ln(p+h) ...)
  Y_zeta(u)  = (theta(x) - x) / sqrt(x)                                 regressor for the zeros of zeta
  Y_chi(u)   = theta(x, chi) / sqrt(x)                                  regressor for the zeros of L(s, chi)
  Y_zeta2, Y_chi2: the same at 2x (theta(2x) - 2x, theta(2x, chi)) / sqrt(x), for the Sophie Germain pattern.
Fit: whiten by differencing (a random-walk noise becomes white), weight by L^(-k/2) (L = ln x, k = number of forms),
FFT with a rectangular window, keep the bins with gamma_lo <= gamma <= gamma_hi, weighted least squares on the
real and imaginary parts, noise variance per bin from the smoothed squared residuals (two passes).
"""
import numpy as np

DU = 1.0 / 2048


def load_cells(paths):
    """Sum the cell arrays of several runs on the global grid index k. Returns (k, X, dict of cell arrays)."""
    zs = [np.load(p) for p in paths]
    kmin = min(int(z["ks"][0]) for z in zs)
    kmax = max(int(z["ks"][-1]) for z in zs)
    k = np.arange(kmin, kmax + 1)
    X = np.zeros(len(k), dtype=np.int64)
    out = {}
    for z in zs:
        o = int(z["ks"][0]) - kmin
        n = len(z["ks"])
        X[o:o + n] = z["X"]
        for key in z.files:
            if key in ("ks", "X", "meta", "next_s"):
                continue
            a = z[key]
            if key not in out:
                out[key] = np.zeros((len(k),) + a.shape[1:])
            out[key][o:o + n] += a
    return k, X, out


def cumulative(cells):
    return {key: np.cumsum(a, axis=0) for key, a in cells.items()}


def char_series(cum, chi_v, q, doubled=False):
    """theta(x, chi) on the grid from the class sums mod q (chi_v: values mod q)."""
    key = "c%d%s" % (q, "_2" if doubled else "")
    return cum[key] @ chi_v


def fit(y, regs, u, j0, j1, glo, ghi, k_forms, nuis=True, smooth_bins=41, return_resid=False):
    """Band-limited WLS of y on the regressors (list of real arrays) over grid indices j0..j1.
    Returns beta, covariance, dict of diagnostics."""
    uu = u[j0:j1 + 1]
    L = uu[:-1]
    w = L ** (-k_forms / 2.0)
    def prep(a):
        a = a[j0:j1 + 1]
        return np.diff(a) * w
    D = prep(y)
    R = [prep(r) for r in regs]
    G = []
    if nuis:
        for g in (uu, uu ** 2, np.exp(-uu / 2), np.exp(-uu)):
            G.append(np.diff(g) * w)
        G.append(w.copy())                     # a constant drift in Y (a linear term in u)
    N = len(D)
    gam = 2 * np.pi * np.fft.rfftfreq(N, d=DU)
    band = (gam >= glo) & (gam <= ghi)
    F = lambda a: np.fft.rfft(a)[band]
    Fd = F(D)
    Xc = np.array([F(r) for r in R] + [F(g) for g in G]).T          # bins x params (complex)
    yr = np.concatenate([Fd.real, Fd.imag])
    Xr = np.concatenate([Xc.real, Xc.imag])
    scale = np.sqrt(np.sum(Xr ** 2, axis=0))
    scale[scale == 0] = 1.0
    Xs = Xr / scale                                            # unit columns (conditioning); undone below
    nb = band.sum()
    var = np.ones(nb)
    for it in range(3):
        W = 1.0 / np.concatenate([var, var])
        XtW = Xs.T * W
        A = XtW @ Xs
        bs = np.linalg.lstsq(A, XtW @ yr, rcond=None)[0]
        res = yr - Xs @ bs
        r2 = res[:nb] ** 2 + res[nb:] ** 2                     # |residual|^2 per bin (2 real dof)
        ker = np.ones(smooth_bins) / smooth_bins
        pad = smooth_bins // 2
        rp = np.concatenate([r2[pad:0:-1], r2, r2[-2:-pad - 2:-1]])
        var = np.convolve(rp, ker, mode="valid") / 2.0
    cov = np.linalg.pinv(A) / np.outer(scale, scale)
    b = bs / scale
    chi2 = float(np.sum(res ** 2 * W))
    dof = len(yr) - len(b)
    diag = dict(nbins=int(nb), chi2_dof=chi2 / dof, gam=gam[band], var=var)
    if return_resid:
        diag["resid"] = res[:nb] + 1j * res[nb:]
        diag["Fd"] = Fd
        diag["Xc"] = Xc
    return b[:len(R)], cov[:len(R), :len(R)], diag
