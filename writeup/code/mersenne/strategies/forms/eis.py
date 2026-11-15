"""Question 4, method A: the eisenstein layout by direct search.

n = a^2 + a b + b^2 with a >= b >= 0, smallest b first (as the page draws it).
Counts the primes drawn and compares the set with {3} u {p = 1 mod 3}.
Checked by eis_check.py (PARI qfbsolve, sympy primes).
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


def main():
    isp = sieve(N)
    primes = np.nonzero(isp)[0]
    bestb = np.full(N + 1, -1, dtype=np.int64)
    besta = np.full(N + 1, -1, dtype=np.int64)
    b = 0
    while 3 * b * b <= N:
        # a >= b and a^2 + a b + b^2 <= N
        amax = (-b + isqrt(b * b - 4 * (b * b - N))) // 2
        while amax * amax + amax * b + b * b > N:
            amax -= 1
        while (amax + 1) ** 2 + (amax + 1) * b + b * b <= N:
            amax += 1
        a = np.arange(b, amax + 1, dtype=np.int64)
        n = a * a + a * b + b * b
        new = bestb[n] < 0
        besta[n[new]] = a[new]
        bestb[n[new]] = b
        b += 1
    drawn = bestb[primes] >= 0
    drawn_set = set(int(p) for p in primes[drawn])
    expected = set(int(p) for p in primes if p == 3 or p % 3 == 1)
    n_drawn_numbers = int((bestb[1:] >= 0).sum())
    res = {
        "pi_N": len(primes),
        "primes_drawn": len(drawn_set),
        "fraction": len(drawn_set) / len(primes),
        "fraction_excluding_3": (len(drawn_set) - 1) / (len(primes) - 1),
        "equals_{3}u{1 mod 3}": drawn_set == expected,
        "primes_1_mod_3": int(sum(1 for p in primes if p % 3 == 1)),
        "primes_2_mod_3": int(sum(1 for p in primes if p % 3 == 2)),
        "numbers_drawn_1..N": n_drawn_numbers,
        "rep_of_3": [int(besta[3]), int(bestb[3])],
        "rep_of_7": [int(besta[7]), int(bestb[7])],
    }
    print(res)
    json.dump(res, open("eis.json", "w"), indent=1)
    np.save("eis_reps.npy", np.stack([besta, bestb]))


if __name__ == "__main__":
    main()
