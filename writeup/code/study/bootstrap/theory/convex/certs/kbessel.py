"""Rigorous K0(x), K1(x) for real x > 0 in arb ball arithmetic, independent of arb's own bessel_k
(python-flint 0.9.0 arb.bessel_k returned a valid but useless ball, [+/- 5.2e-28] for K0(x) = 3.7e-68 at
x = 6 pi e^{67/32} ~ 153.04 given as a 320-bit ball, see arbrepro.py; here we get full precision).

x >= 120 : asymptotic expansion (L = min(2x, 400) terms, relative error <~ e^{-2x})  K_nu(x) = sqrt(pi/(2x)) e^{-x} [ sum_{k<L} a_k(nu)/x^k + R_L ],
           a_k(nu) = prod_{i=1..k} (4nu^2-(2i-1)^2) / (k! 8^k),
           |R_L| <= |a_L(nu)|/x^L  (DLMF 10.40.10: real nu, x > 0, L >= |nu| - 1/2).
x <  120 : power series (DLMF 10.31.1, 10.25.2), all terms in arb at raised precision,
           K0 = -ln(x/2) I0 + sum_k psi(k+1) q^k/(k!)^2,
           K1 = 1/x + ln(x/2) I1 - (x/4) sum_k (psi(k+1)+psi(k+2)) q^k/(k!(k+1)!),   q = x^2/4,
           summed to k = K with K+1 >= x (then consecutive-term ratios are <= 0.35 for each of the series,
           and the remainder of each positive series is <= its last term).
NOTE: in the series branch the radius of the input ball is amplified by ~e^{2x}; callers must carry
~2.9x more bits than the accuracy they need (certify.py works at 800 bits).
"""
from flint import arb, ctx

_A = {}


def a_coef(nu, k):
    key = (nu, k)
    if key not in _A:
        num = 1
        for i in range(1, k + 1):
            num *= (4 * nu * nu - (2 * i - 1) ** 2)
        den = 1
        for i in range(1, k + 1):
            den *= i * 8
        _A[key] = (num, den)
    return _A[key]


def K01(x, L=None):
    """x: positive arb (may have a radius; then the result encloses K over the ball). returns (K0, K1)."""
    assert x.lower() > 0
    if x.lower() >= 120:
        if L is None:                     # near-optimal truncation: terms decrease up to k ~ 2x
            L = min(int(2 * float(x.lower())), 400)
        pre = (arb.pi() / (2 * x)).sqrt() * (-x).exp()
        out = []
        for nu in (0, 1):
            s = arb(0)
            xk = arb(1)
            for k in range(L):
                num, den = a_coef(nu, k)
                s += arb(num) / arb(den) / xk
                xk = xk * x
            num, den = a_coef(nu, L)
            rem = abs(arb(num) / arb(den)) / xk
            s += arb(0, rem.upper())
            out.append(pre * s)
        return out[0], out[1]
    # power series at raised precision
    old = ctx.prec
    try:
        ctx.prec = old + int(3 * float(x.upper())) + 40
        q = x * x / 4
        Kmin = max(int(float(x.upper())) + 2, 6)          # K+1 >= x  (ratio bounds below)
        # stop when the next term is below 2^-(old+20) e^{-2x} (the sums are ~e^{x}, the result ~e^{-x})
        target = arb(2) ** (-(old + 20)) * (-2 * x.upper()).exp()
        I0 = arb(0); I1 = arb(0); S0 = arb(0); S1 = arb(0)
        term = arb(1)          # q^k/(k!)^2
        term1 = arb(1)         # q^k/(k!(k+1)!)
        psi = -arb.const_euler()                 # psi(1)
        k = 0
        while True:
            psi_k2 = psi + arb(1) / (k + 1)      # psi(k+2)
            I0 += term
            I1 += term1
            S0 += psi * term
            S1 += (psi + psi_k2) * term1
            term = term * q / ((k + 1) * (k + 1))
            term1 = term1 * q / ((k + 1) * (k + 2))
            psi = psi_k2
            k += 1
            if k > Kmin and term.upper() < target.lower():
                break
        # remainder bounds (terms from index k on): consecutive ratios <= 0.35 because q/(k+1)^2 <= 1/4
        # (k+1 >= x) and psi(k+2)/psi(k+1) <= 1.37 (k >= 2), all terms positive for k >= 2:
        # remainder <= next/(1-0.35) <= 1.6*next; factor 2 extra safety on the psi-weighted sums.
        nt, nt1 = term, term1
        r_I0 = (nt * arb("1.6")).upper()
        r_I1 = (nt1 * arb("1.6")).upper()
        r_S0 = (nt * abs(psi) * arb("1.6") * arb(2)).upper()
        r_S1 = (nt1 * (abs(psi) * 2 + 1) * arb("1.6") * arb(2)).upper()
        I0 += arb(0, r_I0); I1 += arb(0, r_I1); S0 += arb(0, r_S0); S1 += arb(0, r_S1)
        I1 = I1 * x / 2
        lg = (x / 2).log()
        K0 = -lg * I0 + S0
        K1 = 1 / x + lg * I1 - x / 4 * S1
        return +K0, +K1
    finally:
        ctx.prec = old
