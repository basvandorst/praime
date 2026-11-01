# Puzzle 963 (writeup/963.md): covers of the first X odd primes by arithmetic progressions of primes.
# Checks each cover file: every line is an arithmetic progression of >= 3 terms, all terms are primes among the first X
# odd primes, and together they contain all of them. With --full it re-derives the lower bound for X = 100 without the
# integer program: branch and bound on the set cover (reduced by dominance), each node pruned by an LP dual bound made
# exactly feasible in fractions. No cover by 25 progressions exists, so 26 is the minimum. The sets include the pairs
# {a, b} with 2b - a prime above 547 (progressions continuing past the range), so the bound holds in both readings.
# Run: python3 p963_verify.py          (a few seconds; needs sympy)
#      python3 p963_verify.py --full   (about 10 minutes; needs scipy too)
import os, sys
from sympy import prime, isprime

here = os.path.dirname(os.path.abspath(__file__))

def check(X, path):
    P = [prime(i) for i in range(2, X + 2)]
    S = set(P)
    aps = [list(map(int, l.split())) for l in open(path) if l.strip()]
    for t in aps:
        assert len(t) >= 3 and all(isprime(p) and p in S for p in t)
        d = t[1] - t[0]
        assert d > 0 and all(t[i + 1] - t[i] == d for i in range(len(t) - 1))
    assert set(p for t in aps for p in t) == S
    print('X = %d: %d progressions cover the first %d odd primes (3 .. %d)' % (X, len(aps), X, P[-1]))
    return P, S

check(100, os.path.join(here, '..', 'data', 'p963_cover_100.txt'))
check(1000, os.path.join(here, '..', 'data', 'p963_cover_1000.txt'))

if '--full' in sys.argv:
    from fractions import Fraction as Fr
    from scipy.optimize import linprog
    import numpy as np
    X = 100
    P = [prime(i) for i in range(2, X + 2)]; S = set(P); pos = {p: i for i, p in enumerate(P)}
    raw = set()
    for i, a in enumerate(P):
        for b in P[i + 1:]:
            d = b - a
            if a - d in S: continue
            t = [a, b]
            while t[-1] + d in S: t.append(t[-1] + d)
            if len(t) >= 3 or isprime(2 * b - a):
                raw.add(sum(1 << pos[p] for p in t))
    SETS = sorted(raw)

    def bits(x):
        while x:
            b = x & -x; yield b.bit_length() - 1; x ^= b

    def reduce(U, F):
        forced = 0; nforced = 0
        while True:
            F = sorted({f & U for f in F if f & U}, key=lambda f: -bin(f).count('1'))
            keep = []
            for f in F:
                if not any(f & ~g == 0 for g in keep): keep.append(f)
            F = keep
            cov = {}
            for j, f in enumerate(F):
                for e in bits(f): cov[e] = cov.get(e, 0) | (1 << j)
            if any(e not in cov for e in bits(U)): return U, F, nforced, False      # uncoverable
            single = [e for e in bits(U) if cov[e] & (cov[e] - 1) == 0]
            if single:
                for e in single:
                    if U >> e & 1:
                        j = cov[e].bit_length() - 1; U &= ~F[j]; nforced += 1
                continue
            # element dominance
            drop = 0
            els = sorted(bits(U), key=lambda e: bin(cov[e]).count('1'))
            for a in range(len(els)):
                e = els[a]
                if drop >> e & 1: continue
                for e2 in els[a + 1:]:
                    if not drop >> e2 & 1 and cov[e] & ~cov[e2] == 0: drop |= 1 << e2
            if drop:
                U &= ~drop; continue
            return U, F, nforced, True

    def lp_bound(U, F):
        els = list(bits(U)); idx = {e: i for i, e in enumerate(els)}
        A = np.zeros((len(F), len(els)))
        for j, f in enumerate(F):
            for e in bits(f): A[j, idx[e]] = 1
        res = linprog(-np.ones(len(els)), A_ub=A, b_ub=np.ones(len(F)), bounds=(0, None), method='highs')
        y = [max(Fr(0), Fr(v).limit_denominator(10000)) for v in res.x]
        worst = max(sum(y[idx[e]] for e in bits(f)) for f in F)
        if worst > 1: y = [v / worst for v in y]
        return sum(y)

    nodes = [0]
    def search(U, F, budget):
        nodes[0] += 1
        U, F, nf, ok = reduce(U, F)
        if not ok: return False
        budget -= nf
        if budget < 0: return False
        if U == 0: return True
        if budget == 0: return False
        lb = lp_bound(U, F)
        if lb > budget: return False
        cov = {}
        for j, f in enumerate(F):
            for e in bits(f): cov.setdefault(e, []).append(j)
        e = min(bits(U), key=lambda e: len(cov[e]))
        # branch: the set covering e is f_1, or not f_1 and f_2, ...
        cands = [F[j] for j in cov[e]]
        rest = F
        for f in cands:
            if search(U & ~f, [g for g in rest if g != f], budget - 1): return True
            rest = [g for g in rest if g != f]                   # later branches exclude earlier choices
        return False

    X, TARGET = 100, 25
    U0 = (1 << X) - 1
    U, F, nf, ok = reduce(U0, SETS)
    print('kernel: %d elements, %d sets, %d forced, LP bound %.3f' % (bin(U).count('1'), len(F), nf, float(lp_bound(U, F)) + nf), flush=True)
    res = search(U0, SETS, TARGET)
    print('cover of the first %d odd primes with <= %d progressions: %s  (%d nodes)' % (X, TARGET, 'EXISTS' if res else 'NONE', nodes[0]))

    assert not res
print('all checks passed')
