"""Puzzle 944.  P a palprime with n digits p1..pn.  Intertwining the digit D gives
N_D = A + D*B with A = sum p_i 100^(n-i) (P spread out with zeros) and B = 10 (100^(n-1) - 1)/99.
So the puzzle asks for a palprime whose spread A starts a prime arithmetic progression with the
fixed difference B.  n must be odd (even-length palindromes are multiples of 11) and n = 1 (mod 3)
makes 3, 7, 13 divide B.  This scans all palprimes with n digits and reports the longest run
D = 0..k of primes (PRP for the 2n-1 digit terms) and the histogram of k.
Usage: p944_scan.py n [report_k]"""
import sys, time
from gmpy2 import mpz, is_prime
n = int(sys.argv[1]); rep = int(sys.argv[2]) if len(sys.argv) > 2 else 4
assert n % 2 == 1
h = (n + 1) // 2
B = mpz(10) * (mpz(100)**(n - 1) - 1) // 99
hist = {}; best = []; t0 = time.time(); npal = 0
import itertools
for first in (1, 3, 7, 9):
    for rest in itertools.product(range(10), repeat=h - 1):
        half = (first,) + rest
        digs = half + half[-2::-1]
        P = 0
        for d in digs: P = P * 10 + d
        if not is_prime(P): continue
        npal += 1
        A = mpz(0)
        for d in digs: A = A * 100 + d
        k = -1
        while k < 9 and is_prime(A + (k + 1) * B): k += 1
        hist[k] = hist.get(k, 0) + 1
        if k >= rep: best.append((k, P)); print(f"k={k} P={P}", flush=True)
print(f"n={n}: {npal} palprimes, histogram of k (D = 0..k prime): {dict(sorted(hist.items()))}, {time.time()-t0:.0f}s")
print("B =", B, "small factors:", [q for q in (2,3,5,7,11,13,37,101) if B % q == 0])
