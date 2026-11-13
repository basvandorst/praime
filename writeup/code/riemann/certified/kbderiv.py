"""Exact integrals of |f^(j)| for the normalised Kaiser-Bessel window f(u) = (1-u^2)^(m/2) I_m(beta sqrt(1-u^2)), |u| < 1.
f(u) = sum_k c_k (1-u^2)^(m+k), c_k = (beta/2)^(m+2k) / (k! (m+k)!), a polynomial series in u; its j-th derivative is
evaluated in mpmath (60 digits) on a grid, and the integral of |f^(j)| by the trapezoid rule on 4,001 points plus 5 %.
For a window of half-width a: int |K^(j)| dv = a^(1-j) * Ij."""
import mpmath as mp
mp.mp.dps = 60
_cache = {}
def Ij(beta, m, j, npts=4001):
    key = (float(beta), m, j)
    if key in _cache: return _cache[key]
    beta = mp.mpf(beta)
    # polynomial coefficients of f in u (even powers), series truncated when terms < 1e-40 relative
    terms = []
    k = 0
    while True:
        ck = (beta / 2) ** (m + 2 * k) / (mp.factorial(k) * mp.factorial(m + k))
        terms.append((m + k, ck))
        if k > 10 and ck < mp.mpf(10) ** -45 * terms[0][1]: break
        k += 1
    deg = 2 * terms[-1][0]
    coef = [mp.mpf(0)] * (deg + 1)
    for N, ck in terms:
        for i in range(N + 1):
            coef[2 * i] += ck * mp.binomial(N, i) * (-1) ** i
    for _ in range(j):
        coef = [coef[d] * d for d in range(1, len(coef))]
    h = mp.mpf(2) / (npts - 1)
    vals = [abs(mp.polyval(coef[::-1], -1 + i * h)) for i in range(npts)]
    I = h * (mp.fsum(vals) - (vals[0] + vals[-1]) / 2)
    _cache[key] = float(I) * 1.05
    return _cache[key]

if __name__ == '__main__':
    import math
    from scipy.special import iv
    for m in (4, 6, 8):
        for j in (4, 6, 8):
            if j > m: continue
            print(f'beta=12 m={m} j={j}: int|f^(j)| = {Ij(12, m, j):.4e}  (K0 = I_m(12) = {iv(m, 12):.4e})')
