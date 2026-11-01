# Puzzle 1020 (also 899 Q2, 897 Q3): squarefree numbers with exactly three prime factors whose reciprocals sum to 1,
# with as few terms as possible (writeup/1020.md).
# A prime p is balanced when the inverses mod p of the cofactors (product of the other two primes) of the terms
# containing p sum to 0 (mod p). If every prime is balanced the sum is an integer; a sum in (0, 2) is then 1.
#   1. bulk: all squarefree pqr <= X0;
#   2. repair every prime v > 13 from the largest down: toggle one or two terms v*a*b with a, b < v, or (v <= 61)
#      an "upward gadget": a term a*v*u with u > v plus terms f*g*u (f, g in F = {2,...,13}) that keep u balanced;
#      the option with the least (change in terms) - LAM*X0*(change in sum) is taken;
#   3. close F exactly: a dynamic programme over the 30030 residue states of F, using toggles of the 20 triples
#      inside F and gadgets (v; pairs f*g) with sum of +-(fg)^-1 = 0 (mod v).
# Usage: python3 p1020_build.py X0 SEED LAM [OUTFILE]      e.g.  python3 p1020_build.py 24600 1 0.7 out.txt
import sys, random
from itertools import combinations
from fractions import Fraction
from math import prod
import numpy as np
from sympy import primerange, isprime
F = [2, 3, 5, 7, 11, 13]; Fset = set(F); M = prod(F)
UPV = 61; UMAX = 20000
def crt(vec):
    s = 0
    for r, f in zip(vec, F):
        Mf = M // f; s = (s + r * Mf * pow(Mf, -1, f)) % M
    return s
def run(X0, seed, lam_mult=1.0, Vmax=20000):
    rnd = random.Random(seed); lam = lam_mult * X0
    P = list(primerange(2, X0 // 6 + 1))
    T = set(); res = {}; nb = {}
    def tog(t, s):                     # t sorted triple, s = +1 add / -1 remove
        a, b, c = t
        (T.add if s > 0 else T.discard)(t)
        for x, y, z in ((a, b, c), (b, a, c), (c, a, b)):
            res[x] = (res.get(x, 0) + s * pow(y * z % x, -1, x)) % x
            nb.setdefault(x, set())
            (nb[x].add if s > 0 else nb[x].discard)((y, z) if y < z else (z, y))
    for i, a in enumerate(P):
        if a * P[i+1] * P[i+2] > X0: break
        for j in range(i+1, len(P)):
            b = P[j]
            if a * b * P[j+1] > X0: break
            for k in range(j+1, len(P)):
                c = P[k]
                if a * b * c > X0: break
                tog((a, b, c), 1)
    S = sum(1.0/(a*b*c) for a, b, c in T)
    for v in sorted((p for p in P if p > 13), reverse=True):
        if not nb.get(v): continue
        r = res.get(v, 0)
        if r == 0: continue
        small = [u for u in P if u < v]
        opts = []
        rinv = pow(r, -1, v)
        for a in small:
            for eps in (1, -1):
                b = (-eps * rinv * pow(a, -1, v)) % v
                if b <= a or b >= v or not isprime(b): continue
                t = (a, b, v)
                present = t in T
                if (eps == 1 and present) or (eps == -1 and not present): continue
                opts.append((eps, [(t, eps)]))
        if v <= UPV:
            # upward gadgets: toggle (a, v, u) with u > v, and keep u balanced with terms (f, g, u), f < g in F
            pairsF = [(f, g) for i, f in enumerate(F) for g in F[i+1:]]
            for a in [x for x in P if x < v and x not in (v,)][:40]:
                for su in (1, -1):
                    # need su * (a*u)^-1 = -r (mod v)  =>  u = (-r * su)^-1 * a^-1 (mod v)
                    c0 = (pow((-r * su) % v, -1, v) * pow(a, -1, v)) % v
                    u = c0
                    while u <= v: u += v
                    tries = 0
                    while u < UMAX and tries < 60:
                        tries += 1
                        if u not in Fset and u != a and isprime(u):
                            t = tuple(sorted((a, v, u)))
                            if (t in T) == (su < 0):
                                need = (-su * pow(a * v % u, -1, u)) % u    # change at u to cancel
                                pres = nb.get(u, set())
                                found = None
                                for k in (1, 2):
                                    for Q in combinations(pairsF, k):
                                        sg = [(-1 if (f, g) in pres else 1) for f, g in Q]
                                        if sum(x * pow(f*g % u, -1, u) for x, (f, g) in zip(sg, Q)) % u == need:
                                            found = [((f, g, u), x) for x, (f, g) in zip(sg, Q)]; break
                                    if found: break
                                if found:
                                    opts.append((su + sum(x for _, x in found), [(t, su)] + found))
                        u += v
        if not opts:
            # two toggles: first any term, then a completing single
            for (a, b) in [pr for pr in nb[v] if pr[1] < v][:60]:
                t1 = (a, b, v); r2 = (r - pow(a*b % v, -1, v)) % v
                if r2 == 0: opts.append((-1, [(t1, -1)])); continue
                r2inv = pow(r2, -1, v)
                for a2 in small:
                    for eps in (1, -1):
                        b2 = (-eps * r2inv * pow(a2, -1, v)) % v
                        if b2 <= a2 or b2 >= v or not isprime(b2): continue
                        t2 = (a2, b2, v)
                        if t2 == t1: continue
                        present = t2 in T
                        if (eps == 1 and present) or (eps == -1 and not present): continue
                        opts.append((eps - 1, [(t1, -1), (t2, eps)]))
                if len(opts) > 30: break
        if not opts: return None
        def score(o):
            ds = sum(e / (t[0]*t[1]*t[2]) for t, e in o[1])
            return o[0] - lam * ds + 0.6 * rnd.random()
        best = min(opts, key=score)
        for t, e in best[1]:
            tog(t, e); S += e / (t[0]*t[1]*t[2])
    # closing DP over residues at F
    target = crt([(-res.get(f, 0)) % f for f in F])
    fi = {f: i for i, f in enumerate(F)}
    groups = []
    for t in combinations(F, 3):
        s = -1 if t in T else 1
        vec = [0]*len(F)
        for x, y, z in ((t[0], t[1], t[2]), (t[1], t[0], t[2]), (t[2], t[0], t[1])):
            vec[fi[x]] = (s * pow(y*z % x, -1, x)) % x
        groups.append([(crt(vec), s, [(t, s)])])
    pairs = list(combinations(F, 2))
    for v in primerange(17, Vmax):
        opts = []
        for k in (1, 2, 3):
            for Q in combinations(pairs, k):
                sg = [(-1 if (min(f, g), max(f, g), v) in T or (f, g, v) in T else 1) for f, g in Q]
                if sum(s * pow(f*g % v, -1, v) for s, (f, g) in zip(sg, Q)) % v: continue
                vec = [0]*len(F)
                for s, (f, g) in zip(sg, Q):
                    vec[fi[f]] = (vec[fi[f]] + s * pow(v*g % f, -1, f)) % f
                    vec[fi[g]] = (vec[fi[g]] + s * pow(v*f % g, -1, g)) % g
                opts.append((crt(vec), sum(sg), [((f, g, v), s) for s, (f, g) in zip(sg, Q)]))
        if opts: groups.append(opts)
    INF = 10**9
    dp = np.full(M, INF, dtype=np.int64); dp[0] = 0
    choice = []
    for g in groups:
        new = dp.copy(); ch = np.full(M, -1, dtype=np.int32)
        for oi, (e, c, ops) in enumerate(g):
            rolled = np.roll(dp, e)
            cand = np.where(rolled >= INF, INF, rolled + c)
            better = cand < new
            new[better] = cand[better]; ch[better] = oi
        choice.append(ch); dp = new
    if dp[target] >= INF // 2: return ('noclose', len(T), S)
    s = target; picks = []
    for gi in range(len(groups) - 1, -1, -1):
        oi = choice[gi][s]
        if oi >= 0:
            e, c, ops = groups[gi][oi]; picks.append(ops); s = (s - e) % M
    assert s == 0
    for ops in picks:
        for t, e in ops: tog(t, e)
    tot = sum(Fraction(1, a*b*c) for a, b, c in T)
    return (tot, len(T), T, S)
if __name__ == '__main__':
    X0, seed, lm = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    out = run(X0, seed, lm)
    if out is None: print('repair failed')
    elif out[0] == 'noclose': print('no close, terms %d, pre-sum %.4f' % (out[1], out[2]))
    else:
        tot, n, T, S0 = out
        print('sum', tot, 'terms', n, 'pre-sum %.4f' % S0, 'largest', max(a*b*c for a, b, c in T))
        if tot == 1 and len(sys.argv) > 4:
            open(sys.argv[4], 'w').write(' '.join(str(a*b*c) for a, b, c in sorted(T, key=lambda t: t[0]*t[1]*t[2])) + '\n')
