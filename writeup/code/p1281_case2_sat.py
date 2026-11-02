# Puzzle 1281 case 2 (writeup/1281.md): is there a string of L digits (1..9) whose 3-digit blocks are pairwise coprime?
# SAT (CaDiCaL via pysat): x[i][d] = digit i is d (exactly one); y[i][v] = block i (digits i, i+1, i+2) equals v,
# for every zero-free v in 111..999: y <-> the three digit literals (exactly one y per block);
# pairwise coprime <=> for every prime p at most one block divisible by p (sequential-counter at-most-one).
# Usage: python3 p1281_case2_sat.py L [first_digits]      (first_digits: optional prefix to fix, e.g. 9 or 99)
import sys, itertools, time
from sympy import factorint, isprime
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType

L = int(sys.argv[1]); prefix = sys.argv[2] if len(sys.argv) > 2 else ''
blocks = {v: factorint(v) for v in range(111, 1000) if '0' not in str(v)}
primes = sorted({p for f in blocks.values() for p in f})
top = [0]
def new():
    top[0] += 1; return top[0]
x = [[0] + [new() for d in range(1, 10)] for i in range(L)]
cls = []
for i in range(L):
    cls.append(x[i][1:]); cls += [[-a, -b] for a, b in itertools.combinations(x[i][1:], 2)]
for i, c in enumerate(prefix): cls.append([x[i][int(c)]])
y = []
for i in range(L - 2):
    yi = {}
    for v in blocks:
        s = str(v); lits = [x[i][int(s[0])], x[i + 1][int(s[1])], x[i + 2][int(s[2])]]
        t = new(); yi[v] = t
        cls.append([t] + [-l for l in lits])                  # all three digits -> y
        cls += [[-t, l] for l in lits]                          # y -> each digit
    y.append(yi)
    cls.append(list(yi.values()))
for p in primes:
    lits = [y[i][v] for i in range(L - 2) for v, f in blocks.items() if p in f]
    if len(lits) > 1:
        enc = CardEnc.atmost(lits=lits, bound=1, top_id=top[0], encoding=EncType.seqcounter)
        top[0] = max(top[0], enc.nv); cls += enc.clauses
t0 = time.time()
s = Cadical153(bootstrap_with=cls); ok = s.solve()
print('L = %d prefix %r: %s (%d vars, %d clauses, %.0f s)' % (L, prefix, 'SAT' if ok else 'UNSAT', top[0], len(cls), time.time() - t0), flush=True)
if ok:
    m = set(l for l in s.get_model() if l > 0)
    out = ''.join(str(next(d for d in range(1, 10) if x[i][d] in m)) for i in range(L))
    bl = [int(out[i:i + 3]) for i in range(L - 2)]
    from math import gcd
    assert all(gcd(a, b) == 1 for a, b in itertools.combinations(bl, 2)) and '0' not in out
    print(out, 'prime' if isprime(int(out)) else 'composite')
