"""obslib.py: the per-interval observables of |L(1/2 + it)| on the grid t_k = 100 + k/20 (k = 0..18000).

Intervals: I_j = [100 + j*theta, 100 + (j+1)*theta), theta = 2*pi, j = 0..142 (143 intervals, all inside [100, 1000]).
Per interval:
  M      = max of ln|L| on the interval: the grid maximum refined by the parabola through it and its two neighbours
           (vertex kept only inside the interval), or the value at an exact interval end (quadratic interpolation
           through the three nearest grid points) if that is larger.
  Mgrid  = the plain grid maximum (for the discretisation check).
  lnZ_b  = ln( (1/theta) int_I |L|^(2b) dt ), trapezoid on the grid with the two end pieces by linear interpolation
           of |L|^(2b) (exact interval length theta), b in BETAS.
  mean_ln = (1/theta) int_I ln|L| dt (same quadrature, ln|L| clipped below at ln 1e-12; diagnostic only).
"""
import numpy as np

T0, H, NPT = 100.0, 0.05, 18001
THETA = 2 * np.pi
NINT = 143
BETAS = (0.5, 1.0, 1.5, 2.0, 3.0)


def grid():
    return T0 + H * np.arange(NPT)


def interval_weights(a, b):
    """indices and trapezoid weights (sum = b - a) for int_a^b f(t) dt with f on the grid, linear interpolation
    in the two partial cells."""
    ka = (a - T0) / H
    kb = (b - T0) / H
    i0 = int(np.floor(ka)); i1 = int(np.floor(kb))
    idx = np.arange(i0, i1 + 2)
    w = np.zeros(len(idx))
    # cell [k, k+1] contributes the integral of the linear interpolant over [max(ka,k), min(kb,k+1)]
    for c in range(i0, i1 + 1):
        lo = max(ka, c); hi = min(kb, c + 1)
        if hi <= lo:
            continue
        # linear interpolant f = f_c (1 - s) + f_{c+1} s, s = x - c
        s0, s1 = lo - c, hi - c
        wc = (s1 - s0) - (s1 ** 2 - s0 ** 2) / 2
        wc1 = (s1 ** 2 - s0 ** 2) / 2
        w[c - i0] += wc * H
        w[c + 1 - i0] += wc1 * H
    return idx, w


_W = None


def weights():
    global _W
    if _W is None:
        _W = [interval_weights(T0 + j * THETA, T0 + (j + 1) * THETA) for j in range(NINT)]
    return _W


def parabola_max(y, k):
    """vertex of the parabola through (k-1, k, k+1); returns (offset in steps, value)."""
    a, b, c = y[k - 1], y[k], y[k + 1]
    den = a - 2 * b + c
    if den >= 0:
        return 0.0, b
    d = 0.5 * (a - c) / den
    d = max(-1.0, min(1.0, d))
    return d, b - 0.25 * (a - c) * d


def observables(absL, betas=BETAS):
    absL = np.asarray(absL, float)
    lnL = np.log(np.maximum(absL, 1e-12))
    out = {'M': np.zeros(NINT), 'Mgrid': np.zeros(NINT), 'argM': np.zeros(NINT), 'mean_ln': np.zeros(NINT)}
    for b in betas:
        out['lnZ_%g' % b] = np.zeros(NINT)
    for j, (idx, w) in enumerate(weights()):
        a = T0 + j * THETA
        bnd = a + THETA
        idx_c = np.clip(idx, 0, NPT - 1)
        vals = absL[idx_c]
        for b in betas:
            out['lnZ_%g' % b][j] = np.log(np.dot(w, vals ** (2 * b)) / THETA)
        out['mean_ln'][j] = np.dot(w, lnL[idx_c]) / THETA
        # grid points strictly inside [a, bnd]
        k0 = int(np.ceil((a - T0) / H - 1e-9)); k1 = int(np.floor((bnd - T0) / H + 1e-9))
        k1 = min(k1, NPT - 1)
        seg = lnL[k0:k1 + 1]
        km = k0 + int(np.argmax(seg))
        out['Mgrid'][j] = lnL[km]
        if 1 <= km <= NPT - 2:
            d, v = parabola_max(lnL, km)
            tv = T0 + (km + d) * H
            if tv < a or tv > bnd:      # vertex outside the interval: keep the grid value (the end is within H)
                d, v = 0.0, lnL[km]
            out['M'][j] = v
            out['argM'][j] = T0 + (km + d) * H - a
        else:
            out['M'][j] = lnL[km]
            out['argM'][j] = T0 + km * H - a
        # the two exact interval ends (not grid points): quadratic interpolation through the 3 nearest grid points
        for x in (a, bnd):
            u = (x - T0) / H; kc = int(min(max(round(u), 1), NPT - 2)); s = u - kc
            ve = lnL[kc] + 0.5 * s * (lnL[kc + 1] - lnL[kc - 1]) + 0.5 * s * s * (lnL[kc + 1] - 2 * lnL[kc] + lnL[kc - 1])
            if ve > out['M'][j]:
                out['M'][j] = ve; out['argM'][j] = x - a
    return out
