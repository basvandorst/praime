"""Predicted Sym^2 family densities at p ~ N^2 (N prime), y = p/N^2, normalization E[sqrt(p) lambda_f(p^2)].
Harmonic (Petersson, moduli c = N^2 u):  M_h(y) = sum_{s,q} K21(y/(s q^2)) / (sqrt(s) phi(s))      (s, q squarefree, coprime)
Natural (harmonic x L(1,Sym^2 f)):       M_nat(y) = sum_{rho} A(rho) K21(rho y),  rho = r/(s q^2),
   A = (6/pi^2) r^{-3/2} w(s) v(q), w(s) = prod l^{3/2}/((l-1)(l^2-1)), v(q) = prod l^2/(l^2-1)
 = closed form (Poisson in r):  24 d_{k2} y^{3/2} - 12/(pi (k-1) sqrt y) sum_{1<=m<4y} nu(m) sqrt(m(4y-m)) U_{k-2}(1 - m/(2y)),
   nu(m) = C0 prod_{l || m} 1/(1-u_l) prod_{l^2 | m} (l^2+l-1)/((l^2-1)(1-u_l)),  u_l = 1/((l-1)(l^2-1)),  C0 = prod (1-u_l).
K21(z) = 4 pi i^{-k} sqrt(z) cos(4 pi z) J_{k-1}(4 pi z)."""
import numpy as np
from scipy.special import jv, eval_chebyu
from sympy import factorint, primerange

PR = np.array(list(primerange(2, 2_000_000)), dtype=float)
C0 = float(np.exp(np.sum(np.log1p(-1.0 / ((PR - 1) * (PR * PR - 1))))))


def nu_local(l, e):
    """local factor nu_l(e) of nu(m), e = v_l(m) (lines r/(s q^2), q any integer: alpha(l,2j,0) = l^j)."""
    if e == 0:
        return 1.0 - 1.0 / ((l - 1) * (l * l - 1))
    v = l * l / (l * l - 1.0)
    val = 1.0 + 1.0 / (l * l - 1) + v * (1 - 1.0 / l) * sum(l ** -j for j in range(1, e // 2 + 1))
    if e % 2 == 1:
        val -= v * l ** (-(e + 1) // 2 - 1)
    return val


def nu(m):
    v = C0
    for l, e in factorint(m).items():
        v *= nu_local(l, e) / nu_local(l, 0)
    return v


_NU = {}


def nu_c(m):
    if m not in _NU:
        _NU[m] = nu(m)
    return _NU[m]


def K21(z, k):
    z = np.asarray(z, dtype=float)
    return 4 * np.pi * (-1) ** (k // 2) * np.sqrt(z) * np.cos(4 * np.pi * z) * jv(k - 1, 4 * np.pi * z)


def M_nat(y, k, hyperbolic=True):
    """natural-weight density; hyperbolic=False gives the Petersson-line version M^P (no -6 sqrt(y)/(k-1))."""
    y = np.atleast_1d(np.asarray(y, dtype=float))
    out = 24.0 * y ** 1.5 if k == 2 else np.zeros_like(y)
    if hyperbolic:
        out = out - 6.0 * np.sqrt(y) / (k - 1)
    mmax = int(np.floor(4 * y.max())) + 1
    acc = np.zeros_like(y)
    for m in range(1, mmax + 1):
        msk = 4 * y > m
        if not msk.any():
            continue
        yy = y[msk]
        acc[msk] += nu_c(m) * np.sqrt(m * (4 * yy - m)) * eval_chebyu(k - 2, 1 - m / (2 * yy))
    return out - 12.0 / (np.pi * (k - 1)) * acc / np.sqrt(y)


def sqfree_pairs(Dmax):
    """(s, q, w(s), v(q), 1/(sqrt(s) phi(s))) for coprime squarefree s, q with s q^2 <= Dmax."""
    out = []
    for q in range(1, int(Dmax ** 0.5) + 1):
        fq = factorint(q)
        if any(e > 1 for e in fq.values()):
            continue
        vq = np.prod([l * l / (l * l - 1.0) for l in fq]) if fq else 1.0
        for s in range(1, Dmax // (q * q) + 1):
            fs = factorint(s)
            if any(e > 1 for e in fs.values()) or any(l in fq for l in fs):
                continue
            ws = np.prod([l ** 1.5 / ((l - 1) * (l * l - 1.0)) for l in fs]) if fs else 1.0
            ph = np.prod([l - 1.0 for l in fs]) if fs else 1.0
            out.append((s, q, ws, vq, 1.0 / (np.sqrt(s) * ph)))
    return out


def M_harm(y, k, Dmax=20000):
    y = np.atleast_1d(np.asarray(y, dtype=float))
    out = np.zeros_like(y)
    for s, q, ws, vq, a in sqfree_pairs(Dmax):
        out += a * K21(y / (s * q * q), k)
    return out


def M_nat_lines(y, k, Dmax=400, R=4000):
    """Direct line sum (check of the closed form): sum_{s,q} sum_{r coprime to sq} (6/pi^2) r^{-3/2} w v K21(r y/(s q^2))
    with the tail sum_{r>R} r^{-3/2} * (mean value 1 of K21) added."""
    from math import gcd
    y = np.atleast_1d(np.asarray(y, dtype=float))
    out = np.zeros_like(y)
    r = np.arange(1, R + 1, dtype=float)
    for s, q, ws, vq, a in sqfree_pairs(Dmax):
        D = s * q * q
        ok = np.array([gcd(int(x), s * q) == 1 for x in r])
        rr = r[ok]
        K = K21(np.outer(rr, y) / D, k)
        tot = (rr[:, None] ** -1.5 * K).sum(axis=0)
        dens = np.prod([1 - 1.0 / l for l in factorint(s * q)]) if s * q > 1 else 1.0
        tot += dens * 2.0 / np.sqrt(R)
        out += (6 / np.pi ** 2) * ws * vq * tot
    return out


if __name__ == "__main__":
    print("C0 = %.6f" % C0, " nu(1..12):", ["%.4f" % nu(m) for m in range(1, 13)])
    ms = np.arange(1, 200001)
    print("mean of nu(m), m <= 2e5: %.5f" % np.mean([nu_c(int(m)) for m in ms]))
    ys = np.array([0.1, 0.2, 0.3, 0.45, 0.7, 1.1, 1.6, 2.3, 2.9])
    for k in (2, 4):
        a = M_nat(ys, k); b = M_nat_lines(ys, k)
        print("k=%d closed form:" % k, np.round(a, 4))
        print("     line sum   :", np.round(b, 4))
        print("     harmonic   :", np.round(M_harm(ys, k), 4))


def M_parts(y, k):
    """(smooth part, class-number part) of M_nat: smooth = 24 d_{k2} y^{3/2} - 6 sqrt(y)/(k-1)."""
    y = np.atleast_1d(np.asarray(y, dtype=float))
    smooth = (24.0 * y ** 1.5 if k == 2 else np.zeros_like(y)) - 6.0 * np.sqrt(y) / (k - 1)
    return smooth, M_nat(y, k) - smooth


def M_nat_nu1(y, k):
    """comparator: M_nat with nu(m) replaced by 1 (no local factors)."""
    y = np.atleast_1d(np.asarray(y, dtype=float))
    out = (24.0 * y ** 1.5 if k == 2 else np.zeros_like(y)) - 6.0 * np.sqrt(y) / (k - 1)
    acc = np.zeros_like(y)
    for m in range(1, int(np.floor(4 * y.max())) + 2):
        msk = 4 * y > m
        if msk.any():
            yy = y[msk]
            acc[msk] += np.sqrt(m * (4 * yy - m)) * eval_chebyu(k - 2, 1 - m / (2 * yy))
    return out - 12.0 / (np.pi * (k - 1)) * acc / np.sqrt(y)


def M_zub(y, k):
    """Zubrilina's weight-k density (Theorem 1 of 2310.07681, squarefree levels, root-number weighted), evaluated at y:
    used only as a GL(2)-shape comparator."""
    y = np.atleast_1d(np.asarray(y, dtype=float))
    P = PR
    alpha = 2 * np.pi * np.exp(np.sum(np.log((1 - P - 2 * P ** 2 + P ** 4) / (P ** 4 - 2 * P ** 2 + P))))
    beta = 2 * np.pi * np.exp(np.sum(np.log((-1 + P ** 2 + P ** 3) / (P * (-1 + P + P ** 2)))))
    gamma = 12 * np.exp(np.sum(np.log(P * (1 + P) / (-1 + P + P ** 2))))
    out = beta / (k - 1) * np.sqrt(y) - (gamma * y if k == 2 else 0)
    for r in range(1, int(2 * np.sqrt(y.max())) + 2):
        f = factorint(r)
        nur = np.prod([1 + l * l / (l ** 4 - 2 * l * l - l + 1.0) for l in f]) if f else 1.0
        msk = 4 * y > r * r
        if msk.any():
            yy = y[msk]
            out[msk] += alpha * (-1) ** (k // 2 - 1) / (k - 1) * nur * np.sqrt(4 * yy - r * r) * eval_chebyu(k - 2, r / (2 * np.sqrt(yy)))
    return out
