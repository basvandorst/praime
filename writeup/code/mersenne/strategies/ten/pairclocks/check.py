"""Independent check of count.py and predict.py, plain integers and sympy.

1. Clocks by definition: ord_q(2) from sympy.n_order; z(q) as the least divisor n of q - (5/q) with
   F_n = 0 mod q, F_n from 2x2 matrix powers (count.py uses odd parts and fast doubling instead).
2. Pair tallies for all d = 2..30 with p <= 10^6, and for d = 2 with p <= 10^7, against count_X.json.
3. Predictions: a second enumeration mod 2^20 * 5 with no special treatment of the deep class;
   it must agree with predict.py within 2^-17.
"""
import json
import sys
from fractions import Fraction as F
import sympy as sp


def matmul(A, B, p):
    return ((A[0] * B[0] + A[1] * B[2]) % p, (A[0] * B[1] + A[1] * B[3]) % p,
            (A[2] * B[0] + A[3] * B[2]) % p, (A[2] * B[1] + A[3] * B[3]) % p)


def fib(n, p):
    R, A = (1, 0, 0, 1), (1, 1, 1, 0)
    while n:
        if n & 1:
            R = matmul(R, A, p)
        A = matmul(A, A, p)
        n >>= 1
    return R[1]


def z_rank(q):
    n = q - sp.legendre_symbol(5, q)
    for t in sp.divisors(n):
        if fib(t, q) == 0:
            return t


def clocks(q):
    return sp.n_order(2, q) % 2 == 0, z_rank(q) % 2 == 0


def tallies(limit, ds):
    flag = {}
    out = {}
    primes = set(sp.primerange(7, limit + 40))
    for d in ds:
        t = {"N": 0, "H2": 0, "H0": 0, "L2": 0, "L0": 0}
        for p in sorted(primes):
            if p > limit:
                break
            if p + d not in primes:
                continue
            for q in (p, p + d):
                if q not in flag:
                    flag[q] = clocks(q)
            (h1, l1), (h2, l2) = flag[p], flag[p + d]
            t["N"] += 1
            t["H2"] += h1 and h2
            t["H0"] += (not h1) and (not h2)
            t["L2"] += l1 and l2
            t["L0"] += (not l1) and (not l2)
        out[d] = t
    return out


def single_plain(q2, q5, e):
    m8 = q2 % 8
    h = F(1) if m8 in (3, 5) else F(0) if m8 == 7 else 1 - F(1, 2 ** (e - 1))
    if q2 % 4 == 3:
        lu = F(1)
    elif q5 in (2, 3):
        lu = F(0)
    else:
        lu = 1 - F(2, 2 ** e)
    return h, lu


def pred_plain(d, K=20):
    M = 1 << K
    cls5 = [c for c in range(1, 5) if (c + d) % 5]
    s = {"H2": F(0), "L2": F(0)}
    hs = {}
    n = 0
    for a in range(1, M, 2):
        b = (a + d) % M
        ea = ((a - 1) & -(a - 1)).bit_length() - 1 if a != 1 else K
        eb = ((b - 1) & -(b - 1)).bit_length() - 1 if b != 1 else K
        for c in cls5:
            h1, l1 = single_plain(a, c, ea)
            h2, l2 = single_plain(b, (c + d) % 5, eb)
            s["H2"] += h1 * h2
            s["L2"] += l1 * l2
            n += 1
    return {k: v / n for k, v in s.items()}


if __name__ == "__main__":
    js = json.load(open(sys.argv[1]))
    ck = js["checkpoints"]
    ok = True
    for limit, ds in ((10 ** 6, list(range(2, 32, 2))), (10 ** 7, [2])):
        i = ck.index(limit)
        mine = tallies(limit, ds)
        for d in ds:
            theirs = {k: js["tally"][str(d)][k][i] for k in ("N", "H2", "H0", "L2", "L0")}
            same = theirs == mine[d]
            ok &= same
            print("p <= %d, d = %2d: %s %s" % (limit, d, mine[d], "same" if same else "DIFFERENT %s" % theirs))
    sys.path.insert(0, ".")
    from predict import pair
    for d in (2, 4, 6, 8, 10, 12, 16, 20, 30):
        a = pair(d, 12)
        b = pred_plain(d, 20)
        for k in ("H2", "L2"):
            diff = abs(float(a[k] - b[k]))
            ok &= diff < 2 ** -17
            print("prediction d=%2d %s: predict.py %s, plain enumeration mod 2^20 %.9f, diff %.1e" % (
                d, k, a[k], float(b[k]), diff))
    print("ALL CHECKS PASS" if ok else "CHECK FAILED")
