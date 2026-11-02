# Conjecture 36 (writeup/conj-36.md): for an even N, T(N) = number of odd primes p = 3, 5, 7, ... tried until N - p
# is prime; the "quality index" is T(N) / ln^2 N.
# 1. Recomputes T(N) for the 13 numbers in the page's table and for the 69 record-holders of
#    ../data/conj36_records.txt (OEIS A025018/A025019): checks N - p prime, N - q composite for all odd primes q < p,
#    T = pi(p) - 1, and that the least primes increase along the list.
# 2. The model: each trial succeeds with probability 2 C2 S(N) / ln N, S(N) = prod_{q | N, q > 2} (q-1)/(q-2);
#    prints the normalised index T * 2 C2 S(N) / ln^2 N of every record (model: at most about 1, approaching 1),
#    and the extreme-value prediction (ln(X/2) + gamma) / (2 C2 ln X) for the record index at scale X.
# 3. Answers Q2 and Q3 from the records: the first N needing more than 982 trials, and the largest index to 4e18.
# Run: python3 conj36_verify.py   (about 20 seconds; needs sympy)
import os
from math import log
from sympy import isprime, primerange, primepi, factorint

C2 = 0.66016181584686957; GAMMA = 0.5772156649
def trials(N):                      # (number of odd primes tried, least prime p); p = 2 only for N = 4
    t = 0
    for p in primerange(2, N):
        if isprime(N - p): return t, p
        if p > 2: t += 1
def S(N):
    s = 1.0
    for q in factorint(N):
        if q > 2: s *= (q - 1) / (q - 2)
    return s

page = [(419911924, 249), (721013438, 277), (1847133842, 283), (35884080836, 408), (599533546358, 482), (76903574497118, 655),
        (184162477860248, 677), (217361316706568, 692), (389965026819938, 734), (1047610575836828, 838),
        (24925556008175266, 958), (31284177910528922, 982), (43181037765133228, 958)]
for N, T in page:
    t, p = trials(N); assert t + 1 == T, (N, t, T)
print("the page's 13 trial counts recomputed: all correct")
here = os.path.dirname(os.path.abspath(__file__))
rec = [tuple(map(int, l.split())) for l in open(os.path.join(here, '..', 'data', 'conj36_records.txt')) if not l.startswith('#')]
assert len(rec) == 69 and all(rec[i][1] < rec[i + 1][1] for i in range(68))
best = (0, None); first_983 = None
print('%4s %22s %6s %6s %7s %6s %7s' % ('#', 'N', 'p', 'T', 'index', 'S(N)', 'norm'))
for i, (N, p) in enumerate(rec):
    t, q = trials(N); t += 1 if p > 2 else 0; assert q == p and t == int(primepi(p)) - 1, (N, p, q, t)
    idx = t / log(N) ** 2; norm = idx * 2 * C2 * S(N)
    if idx > best[0]: best = (idx, N)
    if t > 982 and first_983 is None: first_983 = (N, p, t)
    if i >= 56: print('%4d %22d %6d %6d %7.4f %6.3f %7.3f' % (i + 1, N, p, t, idx, S(N), norm))
print('all 69 records consistent (N - p prime, no smaller odd prime works, least primes increasing)')
print('Q2: first record needing more than 982 trials: N = %d (least prime %d, %d trials)' % first_983)
print('Q3: largest index among the records to 4e18: %.4f at N = %d' % best)
for X in (1e9, 1e12, 1e15, 4e18, 1e30, 1e100):
    print('model: record index expected at scale %.0e: %.3f' % (X, (log(X / 2) + GAMMA) / (2 * C2 * log(X))), ' limit 1/(2 C2) = %.4f' % (1 / (2 * C2)) if X == 1e100 else '')
print('all checks passed')
