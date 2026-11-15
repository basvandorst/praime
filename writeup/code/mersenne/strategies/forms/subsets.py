"""Question 6 support, method A (direct search): how the pictures nest.

  P=4, 16, 64 inside P=1: for p = 1 mod 4 write p = a^2 + b^2 (b <= a, as
     drawn for P=1); let e be the even one of a, b.  Claim
     p = x^2 + 16 y^2 <=> 4 | e <=> p = 1 mod 8,  p = x^2 + 64 y^2 <=> 8 | e.
  P=12, 27 inside P=3: p = a^2 + 3 b^2 (as drawn for P=3).  Claim
     p = x^2 + 12 y^2 <=> b even,  p = x^2 + 27 y^2 <=> 3 | b.
  P=58: 2x^2 + 29 (the other form class of discriminant -232 at y = 1) is
     prime for x = 0..28; the residue classes mod 232 of x^2 + 58 y^2.
Checked by subsets_check.py (PARI qfbsolve).
"""
import json
from math import isqrt

import numpy as np

N = 500_000


def sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


ISP = sieve(N)
PR = np.nonzero(ISP)[0]


def reps(P):
    besta = np.full(N + 1, -1, dtype=np.int64)
    bestb = np.full(N + 1, -1, dtype=np.int64)
    b = 0
    while P * b * b <= N:
        a = np.arange(isqrt(N - P * b * b) + 1, dtype=np.int64)
        n = a * a + P * b * b
        new = bestb[n] < 0
        besta[n[new]], bestb[n[new]] = a[new], b
        b += 1
    return besta, bestb


def main():
    res = {}
    a1, b1 = reps(1)
    _, b16 = reps(16)
    _, b64 = reps(64)
    q = PR[PR % 4 == 1]
    a, b = a1[q], b1[q]
    e = np.where(a % 2 == 0, a, b)
    r16, r64 = b16[q] >= 0, b64[q] >= 0
    res["P=1 primes 1 mod 4"] = int(len(q))
    res["4|e count"] = int((e % 4 == 0).sum())
    res["P=16 count"] = int(r16.sum())
    res["P=16 <=> 4|e"] = bool(np.array_equal(r16, e % 4 == 0))
    res["P=16 <=> p=1 mod 8"] = bool(np.array_equal(r16, q % 8 == 1))
    res["8|e count"] = int((e % 8 == 0).sum())
    res["P=64 count"] = int(r64.sum())
    res["P=64 <=> 8|e"] = bool(np.array_equal(r64, e % 8 == 0))
    res["P=64 primes on rows b=0 mod 8"] = int(((e % 8 == 0) & (e == b)).sum())
    res["P=64 primes on columns a=0 mod 8"] = int(((e % 8 == 0) & (e == a)).sum())

    a3, b3 = reps(3)
    _, b12 = reps(12)
    _, b27 = reps(27)
    q3 = PR[(PR % 3 == 1)]
    bb = b3[q3]
    assert (bb >= 0).all()
    res["P=3 primes 1 mod 3"] = int(len(q3))
    res["P=12 <=> b even"] = bool(np.array_equal(b12[q3] >= 0, bb % 2 == 0))
    res["P=12 count"] = int((b12[q3] >= 0).sum())
    res["P=27 <=> 3|b"] = bool(np.array_equal(b27[q3] >= 0, bb % 3 == 0))
    res["P=27 count"] = int((b27[q3] >= 0).sum())
    res["P=27 share of P=3 primes (p=1 mod 3)"] = float((bb % 3 == 0).mean())

    run = 0
    while ISP[2 * run * run + 29]:
        run += 1
    res["2x^2+29 prime for x = 0..(run-1); run"] = run
    a58, b58 = reps(58)
    good = PR[(232 % PR) != 0]
    cls = sorted(set(int(p % 232) for p in good if b58[p] >= 0))
    res["P=58 classes mod 232"] = cls
    res["P=58 number of classes"] = len(cls)
    # P = 5: classes mod 20 with (-5/p) = 1 (Euler criterion), classes of
    # x^2 + 5y^2 and of 2x^2 + 2xy + 3y^2 (direct search)
    good5 = [int(p) for p in PR if 20 % p]
    res["P=5 classes with (-5/p)=1"] = sorted({p % 20 for p in good5
                                               if pow((-5) % p, (p - 1) // 2, p) == 1})
    _, b5 = reps(5)
    res["P=5 classes of x^2+5y^2"] = sorted({p % 20 for p in good5 if b5[p] >= 0})
    hit = np.zeros(N + 1, dtype=bool)
    for x in range(0, isqrt(N) + 1):
        for y in range(-isqrt(N), isqrt(N) + 1):
            v = 2 * x * x + 2 * x * y + 3 * y * y
            if v <= N:
                hit[v] = True
    res["P=5 classes of 2x^2+2xy+3y^2"] = sorted({p % 20 for p in good5 if hit[p]})
    res["P=5: every prime with (-5/p)=1 is in exactly one of the two forms"] = all(
        (b5[p] >= 0) != bool(hit[p]) for p in good5 if pow((-5) % p, (p - 1) // 2, p) == 1)
    print(json.dumps(res, indent=1))
    json.dump(res, open("subsets.json", "w"), indent=1)


if __name__ == "__main__":
    main()
