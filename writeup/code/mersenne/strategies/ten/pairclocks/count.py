"""Count prime pairs (p, p+d), 7 <= p <= X, d = 2..30, by the parity of two clocks of each member.

  H(q): ord_q(2) even  <=> 2^m != 1 mod q, m = odd part of q - 1
  L(q): z(q) even      <=> F_m != 0 mod q, m = odd part of q - (5/q)
Segmented numpy sieve; vectorised square-and-multiply and Fibonacci fast doubling.  Below 2^31 a
product of two residues stays below 2^62; above, mulmod splits one factor into 17-bit halves
(valid for q < 2^34, every partial product below 2^52).  Prints a summary and writes count_X.json (tallies at checkpoints up to X; the last one is X).
Usage: python3 count.py [X]
"""
import json
import sys
import time
import numpy as np

X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 9
DS = list(range(2, 32, 2))
PAD = 32
SEG = 10 ** 7
CK = sorted(set(int(round(10 ** (3 + k / 10))) for k in range(0, 71)) | {10 ** j for j in range(3, 11)} | {X})
CK = [c for c in CK if c <= X]
assert CK[-1] == X


def small_primes(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


SP = small_primes(int((X + PAD) ** 0.5) + 2)


def sieve(lo, hi):
    s = np.ones(hi - lo, bool)
    for p in SP:
        if p * p >= hi:
            break
        st = max(p * p, (lo + p - 1) // p * p)
        s[st - lo::p] = False
    if lo <= 1:
        s[:2 - lo] = False
    return s


B17 = np.uint64(17)
M17 = np.uint64((1 << 17) - 1)


def mulmod(a, b, p):
    if X < 2 ** 31:
        return a * b % p
    hi = b >> B17
    lo = b & M17
    return ((a * hi % p) << B17) % p + a * lo % p


def odd_part(n):
    low = n & (~n + np.uint64(1))
    return n // low, np.log2(low.astype(np.float64)).astype(np.int64)


def clock_h(p):
    m, _ = odd_part(p - np.uint64(1))
    res = np.ones_like(p)
    base = np.full_like(p, 2) % p
    e = m.copy()
    while True:
        bit = (e & np.uint64(1)).astype(bool)
        res = np.where(bit, mulmod(res, base, p) % p, res)
        e >>= np.uint64(1)
        if not e.any():
            break
        base = mulmod(base, base, p) % p
    return res != 1


def clock_l(p):
    r5 = p % np.uint64(5)
    split = (r5 == 1) | (r5 == 4)
    n = np.where(split, p - np.uint64(1), p + np.uint64(1))
    m, _ = odd_part(n)
    a = np.zeros_like(p)
    b = np.ones_like(p)
    for k in range(int(m.max()).bit_length() - 1, -1, -1):
        t = (2 * b + p - a) % p
        c = mulmod(a, t, p) % p
        dd = (mulmod(a, a, p) + mulmod(b, b, p)) % p
        bit = ((m >> np.uint64(k)) & np.uint64(1)).astype(bool)
        a, b = np.where(bit, dd, c), np.where(bit, (c + dd) % p, dd)
    return a != 0


def model(q):
    """Single-prime model probabilities (H, L) as floats, from q mod 8, v2(q-1) and q mod 5."""
    q = q.astype(np.int64)
    low = (q - 1) & -(q - 1)
    e = np.log2(low.astype(np.float64))
    m8 = q % 8
    h = np.where((m8 == 3) | (m8 == 5), 1.0, np.where(m8 == 7, 0.0, 1 - 2.0 ** (1 - e)))
    q5 = q % 5
    lu = np.where(q % 4 == 3, 1.0, np.where((q5 == 2) | (q5 == 3), 0.0, 1 - 2.0 ** (1 - e)))
    return h, lu


def main():
    t0 = time.time()
    keys = ["N", "H2", "H0", "L2", "L0", "EH2", "EH0", "EL2", "EL0", "VH2", "VH0", "VL2", "VL0"]
    tally = {d: {k: np.zeros(len(CK), np.float64 if k[0] in "EV" else np.int64) for k in keys} for d in DS}
    nk = 8 * 12 * 5
    cls_n = np.zeros(nk, np.int64)
    cls_h = np.zeros(nk, np.int64)
    cls_l = np.zeros(nk, np.int64)
    ck = np.array(CK, np.int64)
    for lo in range(0, X + 1, SEG):
        hi = min(lo + SEG, X + 1)
        s = sieve(lo, hi + PAD)
        P = (np.nonzero(s)[0] + lo).astype(np.uint64)
        P = P[P >= 7]
        if len(P) == 0:
            continue
        H = clock_h(P)
        L = clock_l(P)
        main_ = P < hi
        Pm = P[main_].astype(np.int64)
        e1 = odd_part(P[main_] - np.uint64(1))[1]
        key = (Pm % 8) * 60 + np.minimum(e1, 11) * 5 + Pm % 5
        cls_n += np.bincount(key, minlength=nk)
        cls_h += np.bincount(key, weights=H[main_], minlength=nk).astype(np.int64)
        cls_l += np.bincount(key, weights=L[main_], minlength=nk).astype(np.int64)
        Pi = P.astype(np.int64)
        for d in DS:
            q = Pm + d
            ok = s[q - lo]
            p1 = Pm[ok]
            j = np.searchsorted(Pi, p1 + d)
            i = np.searchsorted(Pi, p1)
            h1, h2, l1, l2 = H[i], H[j], L[i], L[j]
            pos = np.searchsorted(p1, ck, side="right")
            mh1, ml1 = model(p1)
            mh2, ml2 = model(p1 + d)
            ex = {"H2": mh1 * mh2, "H0": (1 - mh1) * (1 - mh2), "L2": ml1 * ml2, "L0": (1 - ml1) * (1 - ml2)}
            for k, w in ex.items():
                c = np.concatenate([[0.0], np.cumsum(w)])
                tally[d]["E" + k] += c[pos]
                c = np.concatenate([[0.0], np.cumsum(w * (1 - w))])
                tally[d]["V" + k] += c[pos]
            for k, w in (("N", np.ones(len(p1), bool)), ("H2", h1 & h2), ("H0", ~h1 & ~h2),
                         ("L2", l1 & l2), ("L0", ~l1 & ~l2)):
                c = np.concatenate([[0], np.cumsum(w)])
                tally[d][k] += c[pos]
        if lo % (10 ** 8) == 0:
            print("  segment %d done, %.0f s" % (lo, time.time() - t0), file=sys.stderr, flush=True)
    out = {"X": X, "checkpoints": CK, "d": DS,
           "tally": {str(d): {k: tally[d][k].tolist() for k in keys} for d in DS},
           "classes": {"n": cls_n.tolist(), "h": cls_h.tolist(), "l": cls_l.tolist(),
                       "key": "(p mod 8)*60 + min(v2(p-1),11)*5 + p mod 5"}}
    with open("count_%d.json" % X, "w") as f:
        json.dump(out, f)
    print("X = %d, primes 7..X, %.0f s" % (X, time.time() - t0))
    print(" d       pairs   H both   H none   L both   L none")
    for d in DS:
        t = tally[d]
        n = t["N"][-1]
        print("%2d %11d %.6f %.6f %.6f %.6f" % (d, n, t["H2"][-1] / n, t["H0"][-1] / n,
                                               t["L2"][-1] / n, t["L0"][-1] / n))


if __name__ == "__main__":
    main()
