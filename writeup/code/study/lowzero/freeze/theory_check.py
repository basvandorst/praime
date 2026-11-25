"""theory_check.py: independent checks (mpmath quadrature, no numpy code shared with the other scripts) of the closed
forms in derivation.md.  r = p^(-1/2), w = r e^{i psi}, Y = -ln|1 - w|, psi uniform.
 (1) E e^{2 b Y} = E|1-w|^{-2b} = 2F1(b, b; 1; 1/p);   b = 1: 1/(1 - 1/p).
 (2) E (r^2 - r cos psi) |1-w|^{-4} = -r^2/(1 - r^2)^2 = -p/(p-1)^2.
 (3) wave-closure PPP shift (theta -> infinity): Delta M = -(1/2) ln E[(1 + d')|1-w|^{-2}],
     d' = (2 ln p / L)(r^2 - r cos psi)/|1-w|^2  =>  -(1/2) ln[(p/(p-1))(1 - 2 ln p/((p-1) L))].
 (4) its first-order agreement with (1/2) ln R_p, R_p = (1 - 1/p)(1 + 2 ln p/((p-1) L')) (L' = L in the comparison).
 (5) twisted second moment algebra: (1 + 1/p) - (2/p)(L' - ln p)/L' = (1 - 1/p)(1 + 2 ln p/((p-1) L')).
 (6) second-order (small beta) form of the quenched per-prime shift: c_p(b) ~ -2 b^2 [Var Y - Var(interval mean of Y)],
     interval 2 pi, compared with the exact spread value at b = 0.25 (both from quadrature)."""
import mpmath as mp

mp.mp.dps = 30


def E(f):
    return mp.quad(lambda s: f(s), [0, mp.pi, 2 * mp.pi]) / (2 * mp.pi)


def main():
    for p in (2, 3, 5, 7):
        r = mp.mpf(p) ** -0.5
        a2 = lambda s: 1 - 2 * r * mp.cos(s) + r * r          # |1-w|^2
        for b in (0.5, 1, 1.5, 2, 3):
            lhs = E(lambda s: a2(s) ** (-b)); rhs = mp.hyp2f1(b, b, 1, mp.mpf(1) / p)
            assert abs(lhs / rhs - 1) < 1e-20, (p, b)
        lhs = E(lambda s: (r * r - r * mp.cos(s)) / a2(s) ** 2)
        assert abs(lhs + mp.mpf(p) / (p - 1) ** 2) < 1e-20
        for L in (9, 11.89, 14):
            dM = -mp.log(E(lambda s: (1 + 2 * mp.log(p) / L * (r * r - r * mp.cos(s)) / a2(s)) / a2(s))) / 2
            closed = -mp.log(mp.mpf(p) / (p - 1) * (1 - 2 * mp.log(p) / ((p - 1) * L))) / 2
            assert abs(dM - closed) < 1e-20
            Rp = (1 - mp.mpf(1) / p) * (1 + 2 * mp.log(p) / ((p - 1) * L))
            tw = (1 + mp.mpf(1) / p) - (mp.mpf(2) / p) * (L - mp.log(p)) / L
            assert abs(tw - Rp) < 1e-25
            print('p=%d L=%5.2f  wave-PPP Delta M = %.5f   (1/2) ln R_p = %.5f   lead (1/2) ln(1-1/p) = %.5f' % (
                p, L, float(dM), float(mp.log(Rp) / 2), float(mp.log(1 - mp.mpf(1) / p) / 2)))
        # (6) small-beta check on an interval of length 2 pi
        lp = mp.log(p); th = 2 * mp.pi; b = mp.mpf('0.25')
        Y = lambda s: -mp.log(mp.sqrt(a2(s)))
        varY = E(lambda s: Y(s) ** 2)
        mean_int = lambda phi: mp.quad(lambda h: Y(phi + h * lp), [0, th]) / th
        var_mean = E(lambda phi: mean_int(phi) ** 2)
        spread = -E(lambda phi: mp.log(mp.quad(lambda h: mp.exp(2 * b * Y(phi + h * lp)), [0, th]) / th))
        print('p=%d  Var Y = %.5f (Li2(1/p)/2 = %.5f)  Var(interval mean) = %.5f  small-beta c(0.25) = %.5f  '
              'spread c(0.25) = %.5f' % (p, float(varY), float(mp.polylog(2, mp.mpf(1) / p) / 2), float(var_mean),
                                         float(-2 * b * b * (varY - var_mean)), float(spread)))
    print('all identities hold')


if __name__ == '__main__':
    main()
