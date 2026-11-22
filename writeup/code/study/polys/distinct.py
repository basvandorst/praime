"""Distinct primes in the best run of consecutive prime values. For lines 'a k b c d e start run' from qq: the longest run of consecutive y with |F(y)| prime and all values distinct
(sympy-free: gmpy2; the final claims are re-checked with sympy). usage: distinct.py FILE [MIN]"""
import sys, gmpy2
MIN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
rows = []
for line in open(sys.argv[1]):
    p = list(map(int, line.split()))
    if len(p) != 8: continue
    a, k, b, c, d, e = p[:6]
    F = lambda y: (((a * y + b) * y + c) * y + d) * y + e
    best = (0, None); run = []; rs = None
    for y in range(-90, 92):
        v = abs(F(y) // k) if y <= 90 else 0
        if y <= 90: assert F(y) % k == 0
        if v > 1 and gmpy2.is_prime(v):
            if not run: rs = y
            run.append(v)
        else:
            if len(set(run)) > best[0]: best = (len(set(run)), rs, len(run))
            run = []
    if best[0] >= MIN: rows.append((best[0], best[2], a, k, b, c, d, e, best[1]))
for r in sorted(set(rows), reverse=True)[:15]: print(r)
