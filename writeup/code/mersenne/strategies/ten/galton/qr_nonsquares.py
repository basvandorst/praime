"""Square-remainder board, 10 rows (odd primes 3..31), numbers up to 10^6 without the factors 3..31: how many are
nonzero squares mod all ten primes, how many of those are real squares, which non-squares are = 1 mod 8 (Lehmer
pseudosquares), and bin 9; the prime bins for 32 <= p <= 10^6.   usage: qr_nonsquares.py"""
import math
from sympy import primerange
X = 10 ** 6; Q = list(primerange(3, 32))
def qr(n, q): return n % q != 0 and pow(n, (q - 1) // 2, q) == 1
cop = [n for n in range(1, X + 1) if all(n % q for q in Q)]
allqr = [n for n in cop if all(qr(n, q) for q in Q)]
sq = [n for n in allqr if math.isqrt(n) ** 2 == n]
print('numbers up to 10^6 without the factors 3..31:', len(cop), '; coin expectation for bin 10: %.1f, bin 9: %.1f' % (len(cop) / 1024, len(cop) * 10 / 1024))
print('in bin 10 (square mod all ten):', len(allqr), '; real squares', len(sq), '; non-squares', len(allqr) - len(sq), '(%.2f of the coin)' % ((len(allqr) - len(sq)) / (len(cop) / 1024)))
print('non-squares in bin 10 that are 1 mod 8 (Lehmer pseudosquares for the primes up to 31):', [n for n in allqr if math.isqrt(n) ** 2 != n and n % 8 == 1])
print('bin 9 for these numbers:', sum(1 for n in cop if sum(qr(n, q) for q in Q) == 9))
P = list(primerange(32, X + 1)); b = [0] * 11
for p in P: b[sum(qr(p, q) for q in Q)] += 1
print('primes 32..10^6 by bin:', b, '; coin bin 10 %.1f (ratio %.2f)' % (len(P) / 1024, b[10] / (len(P) / 1024)))
