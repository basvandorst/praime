"""Puzzle 944, second stage.  Reads palprimes P (one per line, survivors of p944_sieve.c), computes
A = P spread with zeros, B = 10 (100^(n-1) - 1)/99 and the largest k with A + D B prime for D = 0..k.
Prints every P with k >= report and the histogram of k.   Usage: p944_ap.py file [report]"""
import sys
from gmpy2 import mpz, is_prime
rep = int(sys.argv[2]) if len(sys.argv) > 2 else 5
hist = {}; m = 0
for line in open(sys.argv[1]):
    s = line.strip()
    if not s.isdigit(): continue
    n = len(s); B = mpz(10) * (mpz(100)**(n - 1) - 1) // 99
    A = mpz(int("0".join(s)))
    k = -1
    while k < 9 and is_prime(A + (k + 1) * B): k += 1
    hist[k] = hist.get(k, 0) + 1; m += 1
    if k >= rep: print(f"k={k} P={s} A={A}", flush=True)
print(f"{m} palprimes tested; histogram of k: {dict(sorted(hist.items()))}")
