import json, mpmath as mp
from chars import chars
mp.mp.dps = 30
out = {}
lg = mp.log(2); g = mp.euler


def Lh(s, q, chi):
    return mp.power(q, -s) * mp.fsum(chi[a] * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q) if chi[a] != 0)


for name, q, v, a, real in chars():
    import cmath
    chi = [mp.mpf(0) if z == 0 else mp.expjpi(mp.mpf(2 * (round(cmath.phase(complex(z)) / (2 * cmath.pi) * 2520) % 2520)) / 2520) for z in v]
    if real: chi = [mp.re(z) for z in chi]
    # L(1) and L'(1) from the Laurent expansion of Hurwitz zeta: zeta(s,x) = 1/(s-1) + gamma0(x) - gamma1(x)(s-1) + ...
    S0 = mp.fsum(chi[b] * mp.stieltjes(0, mp.mpf(b) / q) for b in range(1, q) if chi[b] != 0)
    S1 = mp.fsum(chi[b] * mp.stieltjes(1, mp.mpf(b) / q) for b in range(1, q) if chi[b] != 0)
    dl1 = -S1 / S0 - mp.log(q)
    # independent check: symmetric difference of the Hurwitz sum at s = 1 +- h (error O(h^2))
    h = mp.mpf('1e-7')
    Lp, Lm = Lh(1 + h, q, chi), Lh(1 - h, q, chi)
    dl1_chk = ((Lp - Lm) / (2 * h)) / ((Lp + Lm) / 2)
    # functional equation: L'/L(0,chi) (odd) or B0 (even) = -ln(q/pi) + gamma + ln 2 - conj(L'/L(1,chi))
    c0 = -mp.log(mp.mpf(q) / mp.pi) + g + lg - mp.conj(dl1)
    if a == 1:
        # exact Hurwitz values at 0: zeta(0,x) = 1/2 - x, zeta'(0,x) = lnGamma(x) - ln(2pi)/2
        L0 = mp.fsum(chi[b] * (mp.mpf(1) / 2 - mp.mpf(b) / q) for b in range(1, q))
        dL0 = mp.fsum(chi[b] * (mp.loggamma(mp.mpf(b) / q) - mp.log(2 * mp.pi) / 2) for b in range(1, q)) - mp.log(q) * L0
        direct = dL0 / L0
        Ca = dl1 - direct + lg
    else:
        e = mp.mpf('1e-12')
        f = lambda s: (Lh(s + e, q, chi) - Lh(s - e, q, chi)) / (2 * e) / Lh(s, q, chi) - 1 / s
        hh = mp.mpf('1e-5'); direct = (f(hh) + f(-hh)) / 2
        Ca = dl1 - direct - lg
    Cb = 2 * mp.re(dl1) + mp.log(mp.mpf(q) / mp.pi) + mp.digamma((1 + mp.mpf(a)) / 2)
    out[name] = dict(q=q, a=a, real=real, dl1=[float(mp.re(dl1)), float(mp.im(dl1))], dl1_chk_diff=float(abs(dl1 - dl1_chk)),
                     c0=[float(mp.re(c0)), float(mp.im(c0))], c0_direct_diff=float(abs(direct - c0)),
                     Ca=[float(mp.re(Ca)), float(mp.im(Ca))], Cb=float(Cb), Ca_Cb_diff=float(abs(Ca - Cb)))
    print(name, 'C_y1', mp.nstr(Ca, 10), 'C_had', mp.nstr(Cb, 12), '|diff| %.1e' % float(abs(Ca - Cb)),
          'dl1', mp.nstr(dl1, 10), 'chk %.1e' % float(abs(dl1 - dl1_chk)), 'c0', mp.nstr(c0, 10),
          'c0 direct diff %.1e' % float(abs(direct - c0)), flush=True)
json.dump(out, open('consts.json', 'w'), indent=1)
