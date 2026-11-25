"""Task 1(c), asymptotic form: is the ratios-conjecture three-point function (tplib) at fixed unfolded distances that of
CUE(N_eff) with the same N_eff, to order 1/L^2, also with bad primes?  CUE(N) at order 1/N^2 equals a random dilation of
the sine process with variance -1/(6N^2) (ident.py), so the claim is
   R3(x; L, S) = R3_sine(x) - (Lambda_S/L^2) f''(1) - (Q_S/L^3) f'''(1) + O(L^-4),   f(rho) = rho^3 R3_sine(rho x),
(the 1/L^3 term: the dilation's third cumulant -6Q/L^3).  We fit L^2 (R3 - R3_sine) = a + b/L + c/L^2 over L = 30..480."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib')
import numpy as np
import loadlibs
bk = loadlibs.load_bklib(); tp = loadlibs.load_tplib(); G = tp.G
PS = tp.PrimeSums(PB=10 ** 4)
def R3(t, s, L, S):
    u, v, w = np.array([0.0]), np.array([t]), np.array([s])
    d = tp.T3_parts(u, v, w, L, tuple(S), PS=PS); T3 = tp.T3_gue(u, v, w) + d['Q'] + d['P']
    T2 = lambda r: -np.sinc(r) ** 2 + bk.C_of_r(G, L, tuple(S), r=np.array([abs(r)]))
    return float((1 + T2(t) + T2(s) + T2(s - t) + T3)[0])
def R3s(t, s, rho=1.0):
    a, b, c = np.sinc(rho * t), np.sinc(rho * (s - t)), np.sinc(rho * s)
    return rho ** 3 * (1 - a * a - b * b - c * c + 2 * a * b * c)
def fder(t, s, n, h=1e-3):
    f = lambda r: R3s(t, s, r)
    if n == 2: return (f(1 + h) - 2 * f(1) + f(1 - h)) / h ** 2
    return (f(1 + 2 * h) - 2 * f(1 + h) + 2 * f(1 - h) - f(1 - 2 * h)) / (2 * h ** 3)
lam = lambda S: 1.573151 + sum(math.log(p) ** 2 / (p - 1) for p in S)
Q = lambda S: 2.315846 - sum(math.log(p) ** 3 / (p - 1) ** 2 for p in S)
Ls = np.array([30., 60., 120., 240., 480.])
print('config (t, s) | S | fitted a (L^-2 coef) vs -Lambda_S f\'\' | fitted b (L^-3 coef) vs -Q_S f\'\'\'')
for (t, s) in ((0.5, 1.0), (0.4, 1.3), (0.8, 1.6), (0.7, 2.0), (1.0, 2.5)):
    for S in ((), (2,), (3,), (2, 3)):
        y = np.array([(R3(t, s, L, S) - R3s(t, s)) * L ** 2 for L in Ls])
        A = np.vstack([np.ones_like(Ls), 1 / Ls, 1 / Ls ** 2]).T; c = np.linalg.lstsq(A, y, rcond=None)[0]
        print('  (%.1f,%.1f) %-6s  a = %+.5f  pred %+.5f  | b = %+.4f  pred %+.4f' % (t, s, str(S), c[0], -lam(S) * fder(t, s, 2), c[1], -Q(S) * fder(t, s, 3)))
