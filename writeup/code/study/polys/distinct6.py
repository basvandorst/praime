"""Distinct primes in the best run of consecutive prime values |N(y)/K|, N = y^6 + B5 y^5 + b4 y^4 + ... + b0.
lines: K B5 b4 b3 b2 b1 b0 start run. usage: distinct6.py FILE..."""
import sys, gmpy2
rows = set()
for fn in sys.argv[1:]:
    for line in open(fn):
        K, B5, b4, b3, b2, b1, b0 = map(int, line.split()[:7])
        F = lambda y: ((((((y + B5) * y + b4) * y + b3) * y + b2) * y + b1) * y + b0)
        best = (0, 0, None); run = []; rs = None
        for y in range(-100, 102):
            v = abs(F(y) // K) if y <= 100 else 0
            if y <= 100: assert F(y) % K == 0
            if v > 1 and gmpy2.is_prime(v):
                if not run: rs = y
                run.append(v)
            else:
                if len(set(run)) > best[0]: best = (len(set(run)), len(run), rs)
                run = []
        rows.add((best[0], best[1], K, B5, b4, b3, b2, b1, b0, best[2]))
for r in sorted(rows, reverse=True)[:12]: print(r)
