"""Round 2 theory: everything from the Hardy-Littlewood correlation functions, by a cluster expansion.

Notation. lam = ln t. For a finite set B of offsets, HL (conjecture) gives the density of n with n + b prime for
all b in B as S(B)/lam^|B|.  Fix a pattern B (offsets of consecutive primes, 0 = min B, s = max B = span) and
let I = {0 < m < s, m not in B} be the interior.  The chance that B is a run of consecutive primes is

    P(B) = sum_{A subset I} (-1)^|A| S(B u A) / lam^(|B|+|A|)                          (inclusion-exclusion)

Divide by S(B)/lam^|B|: the remaining sum is the avoidance probability of a point process on I with correlation
functions mu(A) = S(B u A)/(S(B) lam^|A|).  Its logarithm is the alternating sum of the factorial cumulants
(Ursell functions) kappa of mu (a standard identity for avoidance probabilities):

    ln [P(B) lam^|B| / S(B)] = - sum_m kappa(m) + sum_{m<m'} kappa(m,m') - sum_{m<m'<m''} kappa(m,m',m'') + ...
    kappa(m)    = rho_m / lam,            rho_m = S(B u {m}) / S(B)
    kappa(m,m') = rho_m rho_m' (R_mm' - 1) / lam^2,   R_mm' = S(B u {m,m'}) S(B) / (S(B u {m}) S(B u {m'}))

Order 0 of this expansion with rho = 1 is the crude Poisson law of round 1; order 1 and order 2 are new here.
No parameter is fitted: every number is a singular series.

Fast singular series: for p > span every residue of the set is distinct, so
    S(H) = prod_{p <= P} (1 - nu_p/p)/(1 - 1/p)^k  *  G_k(P),   G_k(P) = prod_{p > P} (1 - k/p)/(1 - 1/p)^k,
with P >= span; G_k is computed once with the primes to 2e7 and the tail sum_{p>2e7} C(k,2)/p^2.
"""
import math
from functools import lru_cache
import numpy as np
from common import small_primes

_P_ALL = small_primes(20_000_000).astype(np.float64)


@lru_cache(maxsize=None)
def G(k, P):
    p = _P_ALL[_P_ALL > P]
    s = np.sum(np.log1p(-k / p) - k * np.log1p(-1.0 / p))
    X = 2e7
    s += -k * (k - 1) / 2 / (X * math.log(X))
    return float(math.exp(s))


def S_many(sets, P=None):
    """singular series of each row of an int array (n, k); rows must be sets (distinct entries)."""
    sets = np.asarray(sets, dtype=np.int64)
    n, k = sets.shape
    span = int((sets.max(1) - sets.min(1)).max()) if n else 0
    if P is None:
        P = max(span, 2)
    pr = small_primes(P)
    logS = np.zeros(n)
    dead = np.zeros(n, dtype=bool)
    for p in pr:
        p = int(p)
        r = np.sort(sets % p, axis=1)
        nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        f = 1 - nu / p
        dead |= f <= 0
        logS += np.log(np.where(f > 0, f, 1.0)) - k * math.log(1 - 1 / p)
    S = np.exp(logS) * G(k, P)
    S[dead] = 0.0
    return S


def S1(H):
    return float(S_many(np.array([sorted(H)]))[0])


# ---------------- the cluster expansion of one pattern ----------------
def pattern_sums(B, order=2):
    """returns (S(B), A1, A2): ln P(B) = ln S(B) - |B| ln lam - A1/lam + A2/lam^2 (order-2 truncation)."""
    B = sorted(B)
    s = B[-1]
    SB = S1(B)
    if SB == 0:
        return 0.0, 0.0, 0.0
    Bset = set(B)
    I = np.array([m for m in range(1, s) if m not in Bset], dtype=np.int64)
    if len(I) == 0:
        return SB, 0.0, 0.0
    Barr = np.array(B, dtype=np.int64)
    one = np.concatenate([np.tile(Barr, (len(I), 1)), I[:, None]], axis=1)
    rho = S_many(one, P=max(s, 2)) / SB
    A1 = float(rho.sum())
    if order < 2:
        return SB, A1, 0.0
    live = I[rho > 0]
    rl = rho[rho > 0]
    if len(live) < 2:
        return SB, A1, 0.0
    i, j = np.triu_indices(len(live), 1)
    two = np.concatenate([np.tile(Barr, (len(i), 1)), live[i, None], live[j, None]], axis=1)
    S2 = S_many(two, P=max(s, 2)) / SB                 # = rho_m rho_m' R_mm'
    A2 = float(np.sum(S2 - rl[i] * rl[j]))             # sum over pairs of rho rho (R - 1)
    return SB, A1, A2


def window_avg(L, W, fn, n=2000):
    u = np.linspace(math.log(L), math.log(L + W), n + 1)
    w = np.exp(u) / u
    return float(np.trapezoid(w * fn(u), u) / np.trapezoid(w, u))


def share_fn(SB, A1, A2, npr, order):
    """share per prime of the pattern, as a function of lam: S/lam^(npr-1) * exp(-A1/lam + A2/lam^2)."""
    if order == 0:
        return lambda lam: SB / lam ** (npr - 1) * np.exp(-(0.0) / lam)
    if order == 1:
        return lambda lam: SB / lam ** (npr - 1) * np.exp(-A1 / lam)
    return lambda lam: SB / lam ** (npr - 1) * np.exp(-A1 / lam + A2 / lam ** 2)


def pattern_cumulants(B, kmax=4):
    """returns (S(B), [A1, A2, ..., Akmax]) with ln[P(B) lam^|B|/S(B)] = sum_k (-1)^k A_k / lam^k,
    A_k = sum over k-subsets of the live interior of lam^k kappa_k; kappa from the set-partition formula:
      kappa3 = s123 - s12 s3 - s13 s2 - s23 s1 + 2 s1 s2 s3
      kappa4 = s1234 - sum_4 s123 s4 - sum_3 s12 s34 + 2 sum_6 s12 s3 s4 - 6 s1 s2 s3 s4
    with s_X = S(B u X)/S(B)."""
    from itertools import combinations
    B = sorted(B)
    s = B[-1]
    SB = S1(B)
    if SB == 0:
        return 0.0, [0.0] * kmax
    Bset = set(B)
    I = np.array([m for m in range(1, s) if m not in Bset], dtype=np.int64)
    Barr = np.array(B, dtype=np.int64)
    P = max(s, 2)

    def Ssub(cols):
        n = cols.shape[0]
        return S_many(np.concatenate([np.tile(Barr, (n, 1)), cols], axis=1), P=P) / SB

    A = [0.0] * kmax
    if len(I) == 0:
        return SB, A
    s1 = Ssub(I[:, None])
    A[0] = float(s1.sum())
    live = I[s1 > 0]; s1 = s1[s1 > 0]; n = len(live)
    if kmax < 2 or n < 2:
        return SB, A
    i, j = np.triu_indices(n, 1)
    s2 = np.zeros((n, n))
    s2[i, j] = Ssub(np.stack([live[i], live[j]], 1)); s2[j, i] = s2[i, j]
    A[1] = float(np.sum(s2[i, j] - s1[i] * s1[j]))
    if kmax < 3 or n < 3:
        return SB, A
    T = np.array(list(combinations(range(n), 3)), dtype=np.int64)
    a, b, c = T[:, 0], T[:, 1], T[:, 2]
    s3v = Ssub(np.stack([live[a], live[b], live[c]], 1))
    k3 = s3v - s2[a, b] * s1[c] - s2[a, c] * s1[b] - s2[b, c] * s1[a] + 2 * s1[a] * s1[b] * s1[c]
    A[2] = float(k3.sum())
    if kmax < 4 or n < 4:
        return SB, A
    s3 = np.zeros((n, n, n))
    for (x, y, z) in ((a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a)):
        s3[x, y, z] = s3v
    Q = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    a, b, c, d = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
    s4 = Ssub(np.stack([live[a], live[b], live[c], live[d]], 1))
    k4 = (s4
          - s3[a, b, c] * s1[d] - s3[a, b, d] * s1[c] - s3[a, c, d] * s1[b] - s3[b, c, d] * s1[a]
          - s2[a, b] * s2[c, d] - s2[a, c] * s2[b, d] - s2[a, d] * s2[b, c]
          + 2 * (s2[a, b] * s1[c] * s1[d] + s2[a, c] * s1[b] * s1[d] + s2[a, d] * s1[b] * s1[c]
                 + s2[b, c] * s1[a] * s1[d] + s2[b, d] * s1[a] * s1[c] + s2[c, d] * s1[a] * s1[b])
          - 6 * s1[a] * s1[b] * s1[c] * s1[d])
    A[3] = float(k4.sum())
    return SB, A


class Pattern:
    """Bahadur form of the avoidance probability (round 2, second version).

    With p_m = rho_m/lam and Z_m = (p_m - X_m)/(1 - p_m) (mean 0), exactly
        E prod_m (1 - X_m) = prod_m (1 - p_m) * E prod_m (1 + Z_m) = prod (1 - p_m) * sum_A E[prod_A Z].
    The connected parts of E[prod_A Z] are, for |A| = 2, 3, 4,
        c_ab / ((1-p_a)(1-p_b)),   -k3_abc / prod(1-p),   +k4_abcd / prod(1-p)
    with c, k3, k4 the joint cumulants of the X's (from the singular series). Exponentiating them:
        ln P(B) = ln S(B) - |B| ln lam + sum ln(1 - p_m) + sum c~ - sum k3~ + sum k4~        (order 1..4)
    Order 1 is "independent interior numbers with HL densities"; the diagonal X^2 = X is exact here, which the
    factorial-cumulant exponential of the first version got wrong (it gave e^(-x-y-xy) instead of 1-x-y at h=6).
    """

    def __init__(self, B, kmax=4):
        from itertools import combinations
        B = sorted(B)
        self.B = B
        s = B[-1]
        self.SB = S1(B)
        self.kmax = kmax
        self.s1 = np.zeros(0)
        if self.SB == 0:
            return
        Bset = set(B)
        I = np.array([m for m in range(1, s) if m not in Bset], dtype=np.int64)
        Barr = np.array(B, dtype=np.int64)
        P = max(s, 2)
        SB = self.SB

        def Ssub(cols):
            n = cols.shape[0]
            return S_many(np.concatenate([np.tile(Barr, (n, 1)), cols], axis=1), P=P) / SB

        if len(I) == 0:
            return
        s1 = Ssub(I[:, None])
        live = I[s1 > 0]; s1 = s1[s1 > 0]; n = len(live)
        self.s1 = s1
        self.pairs = self.trip = self.quad = None
        if kmax >= 2 and n >= 2:
            i, j = np.triu_indices(n, 1)
            s2v = Ssub(np.stack([live[i], live[j]], 1))
            s2 = np.zeros((n, n)); s2[i, j] = s2v; s2[j, i] = s2v
            self.pairs = (i, j, s2v)
            self.s2 = s2
        if kmax >= 3 and n >= 3:
            T = np.array(list(combinations(range(n), 3)), dtype=np.int64)
            a, b, c = T[:, 0], T[:, 1], T[:, 2]
            s3v = Ssub(np.stack([live[a], live[b], live[c]], 1))
            self.trip = (a, b, c, s3v)
            if kmax >= 4 and n >= 4:
                s3 = np.zeros((n, n, n))
                for (x, y, z) in ((a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a)):
                    s3[x, y, z] = s3v
                Q = np.array(list(combinations(range(n), 4)), dtype=np.int64)
                a, b, c, d = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
                s4v = Ssub(np.stack([live[a], live[b], live[c], live[d]], 1))
                self.quad = (a, b, c, d, s4v, s3)

    def log_avoid(self, lam, order):
        """ln of the avoidance probability of the interior at scale lam (scalar)."""
        if self.SB == 0:
            return -np.inf
        p = self.s1 / lam
        out = float(np.sum(np.log1p(-p)))
        if order >= 2 and self.pairs is not None:
            i, j, s2v = self.pairs
            c = s2v / lam ** 2 - p[i] * p[j]
            out += float(np.sum(c / ((1 - p[i]) * (1 - p[j]))))
        if order >= 3 and self.trip is not None:
            a, b, cc, s3v = self.trip
            m2 = self.s2 / lam ** 2
            k3 = s3v / lam ** 3 - m2[a, b] * p[cc] - m2[a, cc] * p[b] - m2[b, cc] * p[a] + 2 * p[a] * p[b] * p[cc]
            out -= float(np.sum(k3 / ((1 - p[a]) * (1 - p[b]) * (1 - p[cc]))))
        if order >= 4 and self.quad is not None:
            a, b, c, d, s4v, s3 = self.quad
            m2 = self.s2 / lam ** 2
            m3 = lambda x, y, z: s3[x, y, z] / lam ** 3
            k4 = (s4v / lam ** 4
                  - m3(a, b, c) * p[d] - m3(a, b, d) * p[c] - m3(a, c, d) * p[b] - m3(b, c, d) * p[a]
                  - m2[a, b] * m2[c, d] - m2[a, c] * m2[b, d] - m2[a, d] * m2[b, c]
                  + 2 * (m2[a, b] * p[c] * p[d] + m2[a, c] * p[b] * p[d] + m2[a, d] * p[b] * p[c]
                         + m2[b, c] * p[a] * p[d] + m2[b, d] * p[a] * p[c] + m2[c, d] * p[a] * p[b])
                  - 6 * p[a] * p[b] * p[c] * p[d])
            out += float(np.sum(k4 / ((1 - p[a]) * (1 - p[b]) * (1 - p[c]) * (1 - p[d]))))
        return out

    def share(self, L, W, order, ngrid=24):
        """share per prime of this pattern of consecutive primes, averaged over the window (weight 1/ln t)."""
        if self.SB == 0:
            return 0.0
        k = len(self.B)
        ug = np.linspace(math.log(L), math.log(L + W), ngrid)
        la = np.array([self.log_avoid(u, order) for u in ug])
        return window_avg(L, W, lambda u: self.SB / u ** (k - 1) * np.exp(np.interp(u, ug, la)))


def share_k(SB, A, npr, order):
    """share per prime as a function of lam, truncated after order `order` of the cumulant series."""
    def f(lam):
        e = 0.0
        for k in range(order):
            e = e + (-1) ** (k + 1) * A[k] / lam ** (k + 1)
        return SB / lam ** (npr - 1) * np.exp(e)
    return f


# ---------------- exact HL pair sum for the variance of counts in [n, n+H) ----------------
def pair_S(h):
    """S({0,h}) = 2 C2 prod_{p | h, p > 2} (p-1)/(p-2) for even h, 0 for odd h."""
    h = np.asarray(h, dtype=np.int64)
    C2x2 = 2 * 0.66016181584686957393
    out = np.where(h % 2 == 0, C2x2, 0.0)
    for p in small_primes(int(h.max()) + 1)[1:]:
        p = int(p)
        out = np.where(h % p == 0, out * (p - 1) / (p - 2), out)
    return out


@lru_cache(maxsize=None)
def pair_sum(H):
    h = np.arange(1, H, dtype=np.int64)
    return float(np.sum((H - h) * (pair_S(h) - 1.0)))


C2 = 0.66016181584686957393


@lru_cache(maxsize=None)
def aligned_sum(H, r0):
    """sum over offset pairs a < b in one interval [n, n+H), n = r0 (mod 10), of c_H * R(b - a) * u(a) u(b):
    u(a) = 1 if n + a is prime to 10; c_H = 1/((1-1/2)(1-1/5))^2; R(h) = HL pair factor without p = 2, 5:
    R(h) = R0 prod_{p | h, p != 2, 5} (p-1)/(p-2), R0 = C2 / ((1-2/5)/(1-1/5)^2)."""
    a = np.arange(H)
    u = (np.gcd(r0 + a, 10) == 1).astype(np.int64)
    h = np.arange(1, H, dtype=np.int64)
    # N_h = sum_a u(a) u(a+h)
    cu = np.correlate(u, u, mode='full')[H:]        # lags 1..H-1
    R = np.full(H - 1, C2 / (0.6 / 0.64))
    for p in small_primes(H)[1:]:
        p = int(p)
        if p == 5:
            continue
        R = np.where(h % p == 0, R * (p - 1) / (p - 2), R)
    return float(np.sum(cu * R) / ((0.5 * 0.8) ** 2))


def var_ratio_aligned(H, L, W):
    """Var/mean of the prime count in the intervals [L + jH, L + (j+1)H) (grid fixed by L, H divisible by 10):
    1 + (2 * aligned_sum - H^2)/(H lam). The one-point part uses p_a = u(a)/(lam (1-1/2)(1-1/5))."""
    assert H % 10 == 0 and L % 10 == 0
    s = aligned_sum(H, 0)
    # diagonal: sum_a p_a (1 - p_a) = H/lam - sum_a p_a^2, sum p_a^2 = (#u) * 6.25 / lam^2 = H*2.5/lam^2;
    # E c^2 - (E c)^2 = H/lam + 2 s/lam^2 - H^2/lam^2 ; ratio = 1 + (2 s - H^2)/(H lam)
    return window_avg(L, W, lambda lam: 1 + (2 * s - H * H) / (H * lam))


def var_ratio_exact(H, L, W):
    """Var/mean = 1 - 1/lam + (2/(H lam)) sum_{h<H} (H-h)(S(0,h)-1)   (HL pairs; no asymptotics)."""
    ps = pair_sum(H)
    return window_avg(L, W, lambda lam: 1 - 1 / lam + 2 * ps / (H * lam))
