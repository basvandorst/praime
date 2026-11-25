"""Per-L-function predictions of Var(s) and <r~> for a real primitive character chi_D, zeros in [t1, t2].
Variants (all parameter-free):
  law      : level(L) of the exact two-point function with no bad prime (second-order deformation formula, base sine)
             + sum over p | D of the wave-closure change (exact wave of p, sine base; p <= 47 exact, p > 47 second order)
  law2     : the second-order formula with the exact two-point function of the full bad set (non-additive), base sine
  slow     : level(L) + sum_p of the slow-wave (dilation) change: Var from E[1/rho], <r~> unchanged
  bblm     : V0 - (2 Lambda_D/L^2)(1+V0) + 6 Q_D (1+V0)/L^3 for Var (BBLM with the bad-prime constants); <r~> = level(L) (no change)
  ignored  : level(L) only (bad primes ignored)
Height average: 4 Gauss nodes in t, weights proportional to the zero density ln(q t/2pi)."""
import sys, os, math, json
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/lib'); sys.path.insert(0, H)
import numpy as np
from scipy.interpolate import CubicSpline
import closure, gapcalc

TAB = H + '/tables.npz'
LGRID = np.round(np.arange(6.5, 12.76, 0.25), 4)
PSMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
LAM0, Q0 = 1.573151, 2.315846


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return [int(x) for x in np.nonzero(s)[0]]


def berry_v(p, L):
    lp = math.log(p)
    def v(u):
        u = np.asarray(u, float); s = 0; m = 1
        while p ** m < 1e16:
            s = s + np.sin(math.pi * m * u * lp / L) ** 2 / (m * m * p ** m); m += 1
        return (2 / math.pi ** 2) * s
    return v


def build(pmax=1400):
    if os.path.exists(TAB):
        return dict(np.load(TAB))
    import level
    B = gapcalc.Base(); m1, V0, R0, _ = gapcalc.stats_with(B, None)
    lev = np.array([level.level(L) for L in LGRID])                       # (nL, 2): dVar, dr
    P = primes_upto(pmax)
    dv = np.zeros((len(P), len(LGRID))); dr = np.zeros_like(dv); dvs = np.zeros_like(dv)
    ph = 2 * math.pi * np.arange(4096) / 4096
    for i, p in enumerate(P):
        for j, L in enumerate(LGRID):
            w = closure.Wave(p, L)
            if p <= 47:
                _, vw, rw, _ = gapcalc.stats_with(B, w); dv[i, j] = V0 - vw; dr[i, j] = R0 - rw
            else:
                a, b = gapcalc.second_order(B, berry_v(p, L)); dv[i, j] = -a; dr[i, j] = -b
            dvs[i, j] = V0 - ((1 + V0) * np.mean(1 / w.rho(0.0, ph)) - 1)
    out = dict(L=LGRID, P=np.array(P), lev=lev, dv=dv, dr=dr, dvs=dvs, V0=V0, R0=R0)
    np.savez(TAB, **out)
    return out


class Predictor:
    def __init__(self):
        T = build(); self.T = T
        self.V0, self.R0 = float(T['V0']), float(T['R0'])
        self.Lg = T['L']; self.pidx = {int(p): i for i, p in enumerate(T['P'])}
        self.lev = [CubicSpline(self.Lg, T['lev'][:, k]) for k in (0, 1)]
        self.sp = {}

    def _s(self, key, p):
        k = (key, p)
        if k not in self.sp: self.sp[k] = CubicSpline(self.Lg, self.T[key][self.pidx[p]])
        return self.sp[k]

    @staticmethod
    def bad(D):
        q = abs(D); out = []; d = 2
        while d * d <= q:
            if q % d == 0:
                out.append(d)
                while q % d == 0: q //= d
            d += 1
        if q > 1: out.append(q)
        return out

    def at_L(self, L, S, law2=False):
        lv, lr = float(self.lev[0](L)), float(self.lev[1](L))
        res = {'ignored': (self.V0 + lv, self.R0 + lr)}
        dv = sum(float(self._s('dv', p)(L)) for p in S); dr = sum(float(self._s('dr', p)(L)) for p in S)
        res['law'] = (self.V0 + lv + dv, self.R0 + lr + dr)
        dvs = sum(float(self._s('dvs', p)(L)) for p in S)
        res['slow'] = (self.V0 + lv + dvs, self.R0 + lr)
        lam = LAM0 + sum(math.log(p) ** 2 / (p - 1) for p in S); Q = Q0 - sum(math.log(p) ** 3 / (p - 1) ** 2 for p in S)
        res['bblm'] = (self.V0 - (2 * lam / L ** 2) * (1 + self.V0) + 6 * Q * (1 + self.V0) / L ** 3, self.R0 + lr)   # dilation: dVar = (1+V0)(sigma^2 - kappa3), sigma^2 = -2 Lambda/L^2, kappa3 = -6 Q/L^3
        if law2:
            import level
            a, b = level.level(L, tuple(S)); res['law2'] = (self.V0 + a, self.R0 + b)
        return res

    def predict(self, D, t1, t2, law2=False):
        q = abs(D); S = self.bad(D)
        x, w = np.polynomial.legendre.leggauss(4)
        t = (t1 + t2) / 2 + (t2 - t1) / 2 * x; Lt = np.log(q * t / (2 * math.pi)); ww = w * Lt; ww /= ww.sum()
        acc = {}
        for L, wt in zip(Lt, ww):
            r = self.at_L(L, S, law2)
            for k, v in r.items():
                a = acc.setdefault(k, np.zeros(2)); a += wt * np.array(v)
        out = {k: {'var': float(v[0]), 'r': float(v[1])} for k, v in acc.items()}
        out['meanL'] = float((ww * Lt).sum()); out['bad'] = S
        return out


if __name__ == '__main__':
    import time; t0 = time.time(); build(); print('tables built in %.0f s' % (time.time() - t0))
    P = Predictor()
    for D in (-4, 5, -163, 1009, -3, 12, 60, -420):
        print(D, P.predict(D, 100, 600))
