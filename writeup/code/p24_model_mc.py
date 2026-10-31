# Puzzle 24: expected number of primes n in [A, B] whose base-b representations (b in BASES), read in base 10, are all prime.
# P(M_b prime | local data) = [M_b coprime to all q <= Q] / (prod_{q<=Q}(1-1/q) * ln M_b)   (Mertens-corrected Cramer).
# Joint over all bases for q | 30030 (exact, per sampled prime n); for 17 <= q <= QMID per base, assumed independent.
# For bases >= 11 a digit >= 10 makes the base-10 reading meaningless: weight 0.
import numpy as np, math, sys
JM = 30030; SMALLQ = [2, 3, 5, 7, 11, 13]
MIDQ = [q for q in range(17, 400) if all(q % d for d in range(2, int(q**.5) + 1))]
cS = np.prod([1 - 1/q for q in SMALLQ]); cM = np.prod([1 - 1/q for q in MIDQ])
def primes_window(lo, L):
    s = np.ones(L, bool); r = int(math.isqrt(lo + L)) + 1
    b = np.ones(r + 1, bool); b[:2] = False
    for i in range(2, int(r**.5) + 1):
        if b[i]: b[i*i::i] = False
    for q in np.nonzero(b)[0]: q = int(q); s[(-lo) % q::q] = False
    return lo + np.nonzero(s)[0].astype(np.int64)
def density(lo, bases, L=4 * 10**6, nwin=4, rng=np.random.default_rng(1)):
    """expected solutions per prime near lo (mean over sampled primes) -> returns (mean weight, primes/number)"""
    P = np.concatenate([primes_window(int(lo * (1 + 0.1 * rng.random())), L) for _ in range(nwin)])
    W = np.ones(len(P)); midf = 1.0
    for b in bases:
        if b == 10:
            W *= 1.0  # n itself is prime by sampling
            continue
        n = P.copy(); rj = np.zeros(len(P), np.int64); pj = 1
        rm = np.zeros((len(MIDQ), len(P)), np.int64); pm = np.ones(len(MIDQ), np.int64); Q = np.array(MIDQ, np.int64)[:, None]
        val = np.zeros(len(P)); pw = 1.0; bad = np.zeros(len(P), bool)
        while n.any():
            d = n % b; n //= b
            if b > 10: bad |= d >= 10
            rj = (rj + d * pj) % JM; pj = pj * 10 % JM
            rm = (rm + d[None, :] * pm[:, None]) % Q; pm = pm * 10 % Q[:, 0]
            val += d * pw; pw *= 10.0
        surv = np.gcd(rj, JM) == 1
        W *= np.where(bad, 0.0, surv / (cS * np.log(val)))
        midsurv = (rm != 0).all(0)
        midf *= midsurv[surv & ~bad].mean() / cM if (surv & ~bad).any() else 0
    return W.mean() * midf
def count(A, B, bases, npts=6):
    xs = np.geomspace(A, B, npts); ds = [density(x, bases) for x in xs]
    grid = np.geomspace(A, B, 300)
    f = np.exp(np.interp(np.log(grid), np.log(xs), np.log(ds))) / np.log(grid)
    return np.trapezoid(f, grid)
if __name__ == "__main__":
    B9 = list(range(2, 10))
    for X in (5.0006e10, 2.8e11, 3.1e12):
        print(f"bases 2..9 (and 10): expected solutions below {X:.3g}: {count(1e8, X, B9):.2f}", flush=True)
