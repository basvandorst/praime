"""Check of claim (A) for Haar unitary matrices (own code, nothing imported from the repo).

Flow: Q(z) = det(1 - zU) = sum_j (-1)^j e_j z^j, coefficient j times w_j = exp(-t j (N-j)).
mu_i = 1/roots of the flowed Q; p_k = sum mu_i^k; target E|p_k|^2 / N.

(1) exact finite-N formula derived here (Diaconis-Shahshahani moments + Mayer expansion):
    E|p_k|^2 = k^2 e^{-2tNk} sum_{lambda |- k} [prod_i e^{2t m_i^2}/m_i] Phi(lambda)^2 / prod_m c_m!
    with Phi = sum over connected graphs on the parts of prod (e^{2t m_i m_j} - 1); valid k <= N.
(2) Monte Carlo over Haar matrices, roots by numpy.
(3) limit formula K_tau(u) (claim A), terms n = 1, 2, 3 by quadrature.
"""
import math, itertools, sys
import numpy as np


def partitions(k, maxpart=None):
    if maxpart is None:
        maxpart = k
    if k == 0:
        yield []
        return
    for m in range(min(k, maxpart), 0, -1):
        for rest in partitions(k - m, m):
            yield [m] + rest


def ursell(ms, t):
    """connected-graph sum of prod (e^{2t m_i m_j}-1) via W(S) = sum_{B containing s0} U(B) W(S\\B)."""
    n = len(ms)
    full = (1 << n) - 1
    W = {}
    for S in range(1, full + 1):
        idx = [i for i in range(n) if S >> i & 1]
        lw = sum(2 * t * ms[a] * ms[b] for a, b in itertools.combinations(idx, 2))
        W[S] = math.exp(lw)
    W[0] = 1.0
    U = {}
    for S in sorted(range(1, full + 1), key=lambda s: bin(s).count("1")):
        low = S & -S
        rest = S ^ low
        acc = W[S]
        sub = rest
        while True:  # B = low | sub, sub proper subset of rest
            if sub != rest:
                B = low | sub
                acc -= U[B] * W[S ^ B]
            if sub == 0:
                break
            sub = (sub - 1) & rest
        U[S] = acc
    return U[full]


def exact_Ep2(N, t, k):
    tot = 0.0
    for lam in partitions(k):
        one = 1.0
        for m in lam:
            one *= math.exp(2 * t * m * m) / m
        cnt = 1
        for m in set(lam):
            cnt *= math.factorial(lam.count(m))
        tot += one * ursell(lam, t) ** 2 / cnt
    return k * k * math.exp(-2 * t * N * k) * tot


def haar(N, rng):
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / math.sqrt(2)
    Qm, R = np.linalg.qr(Z)
    d = np.diag(R)
    return Qm * (d / abs(d))


def mc_Ep2(N, t, ks, reps, rng):
    j = np.arange(N + 1)
    w = np.exp(-t * j * (N - j))
    acc = np.zeros(len(ks))
    acc2 = np.zeros(len(ks))
    off = 0.0
    for _ in range(reps):
        lam = np.linalg.eigvals(haar(N, rng))
        c = np.poly(lam)  # c[j] = (-1)^j e_j  (coefficient of x^{N-j})
        q = c * w          # flowed coefficients of z^j in Q(z)
        r = np.roots(q[::-1])  # np.roots wants highest power first: z^N coefficient is q[N]
        mu = 1.0 / r
        off = max(off, np.max(np.abs(np.abs(mu) - 1)))
        for a, k in enumerate(ks):
            v = abs(np.sum(mu ** k)) ** 2
            acc[a] += v
            acc2[a] += v * v
    m = acc / reps
    se = np.sqrt(np.maximum(acc2 / reps - m * m, 0) / reps)
    return m, se, off


def K_limit(tau, u, nq=400):
    """n = 1, 2, 3 terms of claim (A)."""
    t1 = u * math.exp(-2 * tau * u * (1 - u))
    # n = 2: (1/2) int_0^u da e^{2tau(a^2+b^2)}/(ab) (e^{2tau ab}-1)^2, b = u - a
    x, wq = np.polynomial.legendre.leggauss(nq)
    a = (x + 1) * u / 2
    wa = wq * u / 2
    b = u - a
    t2 = 0.5 * np.sum(wa * np.exp(2 * tau * (a * a + b * b)) / (a * b) * np.expm1(2 * tau * a * b) ** 2)
    # n = 3: (1/6) int over simplex a+b+c=u; Phi_3 = f_ab f_bc + f_ab f_ac + f_ac f_bc + f_ab f_bc f_ac
    t3 = 0.0
    xs, ws = np.polynomial.legendre.leggauss(120)
    for xa, wa1 in zip(xs, ws):
        aa = (xa + 1) * u / 2
        rem = u - aa
        bb = (xs + 1) * rem / 2
        wb = ws * rem / 2
        cc = rem - bb
        fab = np.expm1(2 * tau * aa * bb)
        fac = np.expm1(2 * tau * aa * cc)
        fbc = np.expm1(2 * tau * bb * cc)
        phi = fab * fbc + fab * fac + fac * fbc + fab * fbc * fac
        g = np.exp(2 * tau * (aa * aa + bb * bb + cc * cc)) / (aa * bb * cc)
        t3 += wa1 * u / 2 * np.sum(wb * g * phi ** 2)
    t3 /= 6
    pre = u * u * math.exp(-2 * tau * u)
    return t1, pre * t2, pre * t3


if __name__ == "__main__":
    rng = np.random.default_rng(20261127)
    N, tau = 10, 3.0
    t = tau / N ** 2
    ks = list(range(1, N + 1))
    ex = [exact_Ep2(N, t, k) for k in ks]
    m, se, off = mc_Ep2(N, t, ks, 40000, rng)
    print(f"N={N} tau={tau}: max ||mu|-1| = {off:.2e}")
    print(" k   exact E|p_k|^2   MC (se)        z-score")
    for k, e, mm, s in zip(ks, ex, m, se):
        print(f"{k:2d}  {e:12.6f}  {mm:10.5f} ({s:.5f})  {(mm - e) / s:+.2f}")
    N, tau = 60, 3.0
    t = tau / N ** 2
    ks = [12, 30, 48]
    m, se, off = mc_Ep2(N, t, ks, 6000, rng)
    print(f"\nN={N} tau={tau}: max ||mu|-1| = {off:.2e}")
    for k, mm, s in zip(ks, m, se):
        u = k / N
        a1, a2, a3 = K_limit(tau, u)
        print(f"u={u:.2f}  MC E|p_k|^2/N = {mm / N:.5f} ({s / N:.5f})   limit n<=1: {a1:.5f}  n<=2: {a1 + a2:.5f}  n<=3: {a1 + a2 + a3:.5f}")
