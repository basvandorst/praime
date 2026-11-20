"""Check of the identity behind the lag of the square arms: for q = 44 and 710 the indicator of 'a is a square mod q'
among the units is (1/c) sum over the real characters (c = number of real characters = number of square roots of a
square), so  mean(non-squares) - mean(squares) = -(c / ((c - 1) phi)) * sum_{real chi != chi0} pi(x, chi)  with
pi(x, chi) = sum_{p <= x} chi(p).  The real characters are built as the characters trivial on the squares.  usage: lag_check.py"""
import math, itertools
import numpy as np
X = 10 ** 7
s = np.ones(X + 1, bool); s[:2] = False
for i in range(2, int(X ** .5) + 1):
    if s[i]: s[i * i::i] = False
P = np.nonzero(s)[0]
for q in (44, 710):
    U = [a for a in range(1, q) if math.gcd(a, q) == 1]; phi = len(U)
    sq = sorted(set(a * a % q for a in U)); c = phi // len(sq)
    cnt = np.bincount(P % q, minlength=q)
    # real characters = functions U -> {+1,-1} multiplicative and trivial on squares: brute force over homomorphisms via generators of U/U^2
    cosets = {}
    for a in U:
        key = min(a * b % q for b in sq); cosets.setdefault(key, a)
    reps = sorted(cosets)                               # c cosets of the squares
    chars = []
    for signs in itertools.product((1, -1), repeat=len(reps)):
        f = {}
        for r, sgn in zip(reps, signs):
            for b in sq: f[r * b % q] = sgn
        if all(f[a * b % q] == f[a] * f[b] for a in U for b in U): chars.append(f)
    sums = [sum(f[a] * int(cnt[a]) for a in U) for f in chars if any(v == -1 for v in f.values())]
    lag_formula = -(c / ((c - 1) * phi)) * sum(sums)
    ns = [cnt[a] for a in U if a not in set(sq)]; ss = [cnt[a] for a in sq]
    print('q = %d: %d real characters (c = %d); sums over p <= 10^7 of the non-principal ones: %s; lag by formula %.4f, measured %.4f' % (q, len(chars), c, sums, lag_formula, np.mean(ns) - np.mean(ss)))
