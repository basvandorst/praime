"""Question 1 and 2, method A (direct search over (a, b)).

For P = 1..100 and every n <= N: record the representation n = a^2 + P b^2
with the smallest b (scan b = 0, 1, 2, ...; first hit wins), exactly as the
page draws it.  Then:
  - count primes drawn (all, and excluding primes dividing 4P),
  - h(-4P) by counting reduced primitive forms,
  - test whether the represented primes (p not dividing 4P) are a union of
    residue classes mod 4P (every class all-in or all-out),
  - compare with Euler's idoneal list.
Writes qf_count.json.  Checked independently by qf_check.py.
"""
import json
from math import gcd, isqrt

import numpy as np

N = 500_000
IDONEAL = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25,
           28, 30, 33, 37, 40, 42, 45, 48, 57, 58, 60, 70, 72, 78, 85, 88, 93}


def sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


def smallest_b_reps(P, n_max=N):
    """besta[n], bestb[n] = representation with smallest b, -1 if none."""
    besta = np.full(n_max + 1, -1, dtype=np.int64)
    bestb = np.full(n_max + 1, -1, dtype=np.int64)
    b = 0
    while P * b * b <= n_max:
        amax = isqrt(n_max - P * b * b)
        a = np.arange(amax + 1, dtype=np.int64)
        n = a * a + P * b * b
        new = bestb[n] < 0
        besta[n[new]] = a[new]
        bestb[n[new]] = b
        b += 1
    return besta, bestb


def class_number_reduced(D):
    """Number of reduced primitive positive definite forms (a,b,c), b^2-4ac=D<0."""
    assert D < 0 and D % 4 in (0, 1)
    h = 0
    forms = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if b < 0 and (a == c):
                continue  # need b >= 0 when a == c (b = -a excluded by range)
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
            forms.append((a, b, c))
        a += 1
    return h, forms


def main():
    isp = sieve(N)
    primes = np.nonzero(isp)[0]
    pi_N = len(primes)
    out = {"N": N, "pi_N": int(pi_N), "rows": []}
    for P in range(1, 101):
        besta, bestb = smallest_b_reps(P)
        rep = bestb[primes] >= 0
        drawn_all = int(rep.sum())
        good = (4 * P) % primes != 0          # p does not divide 4P
        drawn_excl = int((rep & good).sum())
        n_excl = int(good.sum())
        special = [int(p) for p in primes[~good]]
        special_rep = [int(p) for p in primes[~good] if bestb[p] >= 0]
        h, forms = class_number_reduced(-4 * P)
        # residue classes mod 4P
        m = 4 * P
        r = primes[good] % m
        cnt_in = np.bincount(r[rep[good]], minlength=m)
        cnt_out = np.bincount(r[~rep[good]], minlength=m)
        coprime = [x for x in range(m) if gcd(x, m) == 1]
        mixed = [x for x in coprime if cnt_in[x] > 0 and cnt_out[x] > 0]
        empty = [x for x in coprime if cnt_in[x] + cnt_out[x] == 0]
        in_classes = [x for x in coprime if cnt_in[x] > 0 and cnt_out[x] == 0]
        # smallest number of primes in a coprime class (to show the test has teeth)
        min_class = int(min(cnt_in[x] + cnt_out[x] for x in coprime))
        ambiguous = all(b == 0 or b == a or a == c for (a, b, c) in forms)
        row = {
            "P": P, "h": h, "forms": forms,
            "drawn_all": drawn_all, "drawn_excl": drawn_excl, "n_excl": n_excl,
            "special_primes": special, "special_drawn": special_rep,
            "fraction_all": drawn_all / pi_N,
            "fraction_excl": drawn_excl / n_excl,
            "chebotarev": 1 / (2 * h),
            "decided": len(mixed) == 0 and len(empty) == 0,
            "n_mixed_classes": len(mixed), "n_classes": len(coprime),
            "n_empty_classes": len(empty), "min_primes_per_class": min_class,
            "rep_classes": in_classes if len(mixed) == 0 else None,
            "idoneal_list": P in IDONEAL,
            "all_reduced_forms_ambiguous": ambiguous,
        }
        out["rows"].append(row)
        print(f"P={P:3d} h={h:2d} drawn={drawn_all:6d} frac={drawn_all/pi_N:.5f} "
              f"excl={drawn_excl/n_excl:.5f} 1/2h={1/(2*h):.5f} "
              f"decided={row['decided']} mixed={len(mixed)}/{len(coprime)} "
              f"idoneal={P in IDONEAL} amb={ambiguous}", flush=True)
    with open("qf_count.json", "w") as f:
        json.dump(out, f)


if __name__ == "__main__":
    main()
