"""Error budget of the family criterion (q prime, all q-1 characters, x = n, a = c n).

Lambda(n) K(0) = (1/phi(q)) { P(x) - zeta zero sum - E_zeta - sum_{chi != chi0} chibar(n) [ L zero sum + E_kappa ] }
Counting bound for primitive chi of conductor q (Bennett-Martin-O'Bryant-Rechnitzer, Math. Comp. 90 (2021), Thm 1.1):
 N(T,chi) = #{rho: 0<beta<1, |gamma| <= T}, |N - (T/pi) log(qT/(2 pi e)) + chi(-1)/4| <= 0.22737 l + 2 log(1+l) - 0.5,
 l = log(q(T+2)/(2 pi)) > 1.567, T >= 5/7.
GRH for conductor q <= 400,000 up to height max(1e8/q, 3.75e7/q + 200) for odd q (Platt, Math. Comp. 85 (2016)).
"""
import math
import numpy as np
from kbm import majorant, K0, Khat, P, K

def NplusL(t, q):
    l = np.log(q * (t + 2) / (2 * math.pi))
    return t / math.pi * np.log(q * t / (2 * math.pi * math.e)) + 0.25 + 0.22737 * l + 2 * np.log(1 + l) - 0.5

def J_L(T, NT, q, c, beta, m, steps=6000):
    """Upper bound of sum_{|gamma| > T} majZ(c|gamma|) for one L-function with NT zeros in |gamma| <= T,
    where majZ = majorant / (a sqrt(2 pi) beta^m). Upper Riemann-Stieltjes sum on a geometric grid."""
    assert c * T > beta * 1.000001
    nu = m + 0.5; p = nu + 1 / 3
    def majZ(t):
        s = c * np.asarray(t, dtype=float); y = np.sqrt(np.clip(s * s - beta * beta, 1e-300, None))
        return np.minimum(1.0 / (2 ** nu * math.gamma(nu + 1)), 0.7858 * y ** (-p))
    tK = T * 1e6
    ts = T * np.exp(np.linspace(0, math.log(tK / T), steps + 1))
    f = np.minimum.accumulate(majZ(ts))
    Np = NplusL(ts, q)
    s = float(np.sum(np.maximum(Np[1:] - NT, 0) * (f[:-1] - f[1:])))
    C = f[-1] * tK ** p
    integral = C * tK ** (1 - p) * ((math.log(q * tK) / math.pi) / (p - 1) + 1 / (math.pi * (p - 1) ** 2) + 1 / p)
    return s + NplusL(tK, q) * f[-1] + integral

def H_platt(q):
    return max(1e8 / q, (7.5e7 if q % 2 == 0 else 3.75e7) / q + 200)

_IJ = {}
def offline_L(x, a, beta, m, j, q):
    """sum over the zeros of one L(s, chi) with |gamma| > H (anywhere in the strip) of |phi~(rho)|."""
    key = (a / x, beta, m, j)
    if key not in _IJ:
        v = np.linspace(-a, a, 400001); hh = v[1] - v[0]; k = K(v, a, beta, m)
        for _ in range(j): k = np.gradient(k, hh)
        _IJ[key] = float(np.sum(np.abs(k))) * hh * 1.1 / a ** (1 - j)   # scales like a^(1-j)
    I = _IJ[key] * a ** (1 - j)
    H = H_platt(q)
    S = (j / math.pi) * H ** (1 - j) * (math.log(q * H) / (j - 1) + 1 / (j - 1) ** 2)
    return math.exp(a / (2 * x)) * I * x ** j * S

if __name__ == '__main__':
    q = 20011; c = 1.09; m = 4; cc = 1.2
    for beta in (9, 10, 11, 12, 12.7, 14):
        T = cc * beta / c
        NT = T / math.pi * math.log(q * T / (2 * math.pi * math.e)) - 8   # pessimistic found count
        J = J_L(T, NT, q, c, beta, m)
        k0 = K0(beta, m)
        for n in (10, 1000, 9999):
            a = c * n
            tail = c * math.sqrt(n) * math.sqrt(2 * math.pi) * beta ** m * J
            off = offline_L(float(n), a, beta, m, 4, q)
            print(f'beta={beta} T={T:.2f} n={n}: tail/K0 = {tail/k0:.3e}  off/K0 = {off/k0:.1e}  half ln n = {0.5*math.log(n):.2f}')
