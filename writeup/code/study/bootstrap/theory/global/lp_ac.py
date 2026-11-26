# LP: absolutely continuous Lambda-relaxed system.  dLambda = e^u (1 + d(u)) du on [u0, U],
# 0 on [0,u0), e^u du on (U, inf).  d piecewise constant on cells, d >= -1, d <= dmax.
# Re H(1/2+it) = b(t) - F(t);  maximise m = min_t (b - F).
import numpy as np, mpmath as mp, sys
from scipy.optimize import linprog
def bfun(t):
    return float(0.5*(mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*t)) - mp.log(mp.pi)) + 0.5/(0.25 + t*t))
def cellint(t, a, c):   # int_a^c cos(tu) e^{u/2} du
    z = 0.5 + 1j*t
    return np.real((np.exp(z*c) - np.exp(z*a))/z)
def run(x0, U, ncell, dmax, Tmax, dt):
    u0 = np.log(x0)
    edges = np.linspace(u0, U, ncell+1)
    ts = np.arange(0, Tmax, dt)
    bs = np.array([bfun(t) for t in ts])
    base = -cellint(ts, 0.0, u0)           # F contribution of the gap [0,u0)
    M = np.stack([cellint(ts, edges[j], edges[j+1]) for j in range(ncell)], axis=1)
    # variables: d_0..d_{n-1}, m ; maximise m  <=>  minimise -m
    # constraint: b - base - M d >= m  <=>  M d + m <= b - base
    A = np.hstack([M, np.ones((len(ts),1))])
    rhs = bs - base
    c = np.zeros(ncell+1); c[-1] = -1
    bounds = [(-1, dmax)]*ncell + [(None, None)]
    r = linprog(c, A_ub=A, b_ub=rhs, bounds=bounds, method='highs')
    return r, edges, ts, bs, base, M
if __name__ == "__main__":
    for x0 in [1.5, 2.0, 2.5]:
        for U in [np.log(20), np.log(100)]:
            r, edges, ts, bs, base, M = run(x0, U, 60, 20.0, 80.0, 0.02)
            print("x0=%.1f U=e^%.2f status=%s  max min(b-F) = %.5f" % (x0, U, r.status, -r.fun if r.status==0 else float('nan')))
            sys.stdout.flush()
