"""Non-additivity for two bad primes.
(1) Wave closure with the additive warp (density 1 + eps_p + eps_q, the explicit formula has no eps_p*eps_q term):
    contrast(S={p,q}) - contrast(p) - contrast(q), against +2 m_p m_q.
(2) The data: interaction coefficient I(p|D) I(q|D) for (2,3), (2,5), (3,5), (2,7) added to the registered design,
    fresh sets and round 142's sets; against the frozen law's pair term and the closure's."""
import os, numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myfit import load, design, ols, PR
from altlib import clo, factor


def wave(p, th, L, chi=1, K=60):
    k = np.arange(1, K + 1)[:, None]
    c = (chi ** k) * p ** (-k / 2.0) * np.cos(k * th[None, :])
    return (c / k).sum(0), -(2 * np.log(p) / L) * c.sum(0)


def clo2(p, q, L, n=1024):
    th = 2 * np.pi * np.arange(n) / n
    Ap, ep = wave(p, th, L); Aq, eq = wave(q, th, L)
    A = Ap[:, None] + Aq[None, :]; e = ep[:, None] + eq[None, :]
    w = 1 + e; X = A + np.log(w)
    m = np.mean(w * X)
    return -(np.mean(w * X * X) - m * m), -m


if __name__ == '__main__':
    print('closure pair term (additive warp):')
    for L in (10.0, 10.7):
        for p, q in [(2, 3), (2, 5), (3, 5), (2, 7)]:
            v2, m2 = clo2(p, q, L)
            vp, mp = clo(p, int(L * 1000)); vq, mq = clo(q, int(L * 1000))
            print('  L=%.1f (%d,%d): pair term %+.5f   2*dmu_p*dmu_q %+.5f   (mean non-additivity %+.1e)' % (
                L, p, q, v2 - vp - vq, 2 * mp * mq, m2 - mp - mq))
    # law's pair term as frozen: -2 dmu_p dmu_q + NA (NA(2,3) ~ +0.0005 at L = 10 per derivation/oned2_out)
    rows = [r for r in load() if r['D'] not in (197, -395)]
    for nm, rr, keyp in [('fresh', rows, 'cum2_3pt')]:
        v = np.array([np.var(r['la'] - np.log(r['Lj']), ddof=1) for r in rr])
        pr = np.array([r['pred'][keyp] for r in rr])
        for (p, q) in [(2, 3), (2, 5), (3, 5)]:
            inter = np.array([1.0 if (r['q'] % p == 0 and r['q'] % q == 0) else 0.0 for r in rr])
            X, nn = design(rr, extra=None)
            Xi = np.column_stack([X, inter - inter.mean()])
            b = ols(Xi, v)[-1]; bp = ols(Xi, pr)[-1]
            rng = np.random.default_rng(3); B = []
            for _ in range(1000):
                ix = rng.integers(0, len(rr), len(rr)); Xb = Xi[ix].copy(); Xb[:, 1:] -= Xb[:, 1:].mean(0)
                B.append(ols(Xb, v[ix])[-1])
            print('%s interaction (%d,%d): n=%d  data %+.4f +- %.4f   frozen law %+.4f' % (nm, p, q, int(inter.sum()), b, np.std(B, ddof=1), bp))
