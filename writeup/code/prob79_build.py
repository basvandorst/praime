# Problem 79: numbers with exactly three prime factors (counted with multiplicity) whose reciprocals sum to 1,
# with as few terms as possible (writeup/prob-79.md).
# For a prime p let R_p = p^2 * (sum of 1/t over the terms t divisible by p), taken mod p^2: a term p^2*q adds q^-1
# and a term p*m (p not dividing m) adds p*m^-1. The sum is an integer iff R_p = 0 (mod p^2) for every p; a cube
# p^3 can never occur (1/p^3 alone has p-valuation -3).
#   1. bulk: all such numbers <= X0;
#   2. repair every prime v > 7 from the largest down: first the residue mod v with terms v^2*a (a < v),
#      then mod v^2 with one or two terms v*a*b (a <= b < v); score (change in terms) - LAM*X0*(change in sum);
#   3. close 2, 3, 5, 7 exactly: a dynamic programme over the 4*9*25*49 = 44100 states of (R_2, R_3, R_5, R_7),
#      using the 16 terms built from 2, 3, 5, 7 and gadgets v*m (m = product of two of 2, 3, 5, 7) that keep v
#      balanced; its cost also counts the change in sum, so that it lands on 1.
# Usage: python3 prob79_build.py X0 SEED LAM [OUTFILE]      e.g.  python3 prob79_build.py 850 7 1.0 out.txt
import sys, random
from itertools import combinations, combinations_with_replacement
from fractions import Fraction
from math import prod
import numpy as np
from sympy import primerange, isprime, factorint
F = [2, 3, 5, 7]; Fset = set(F); MODS = [p*p for p in F]; M = prod(MODS)
def crt(vec):
    s = 0
    for r, m in zip(vec, MODS):
        Mm = M // m; s = (s + r * Mm * pow(Mm, -1, m)) % M
    return s
def contrib(p, t):
    """contribution of term t (int) to R_p (mod p^2); p | t."""
    if t % (p*p) == 0: return pow(t // (p*p), -1, p*p)
    m = t // p
    return (p * pow(m % p, -1, p)) % (p*p)
def run(X0, seed, lam_mult=0.7, Vmax=6000):
    rnd = random.Random(seed); lam = lam_mult * X0
    T = set(); R = {}
    def tog(t, s):
        (T.add if s > 0 else T.discard)(t)
        for p in factorint(t):
            R[p] = (R.get(p, 0) + s * contrib(p, t)) % (p*p)
    for t in range(8, X0 + 1):
        f = factorint(t)
        if sum(f.values()) == 3 and len(f) > 1: tog(t, 1)
    S = sum(1.0/t for t in T)
    P = list(primerange(2, X0 // 4 + 1))
    for v in sorted((p for p in P if p > F[-1]), reverse=True):
        r = R.get(v, 0)
        if r == 0: continue
        small = [u for u in P if u < v]
        opts = []
        def level1_options(r1, excl=()):
            # toggles of v*a*b (a <= b < v) changing R_v by +- v*(ab)^-1; need sum = -r1 (r1 multiple of v)
            out = []; y = (r1 // v) % v
            if y == 0: return [(0, [])]
            for a in small:
                for eps in (1, -1):
                    # eps*(ab)^-1 = -y (mod v)  =>  b = -eps * (y*a)^-1 (mod v)
                    b = (-eps * pow(y * a % v, -1, v)) % v
                    if b < a or b >= v or not isprime(b): continue
                    t = v * a * b
                    if t in excl: continue
                    if (eps == 1) == (t in T): continue
                    out.append((eps, [(t, eps)]))
            return out
        if r % v == 0:
            opts = [o for o in level1_options(r) if o[1]]
            if not opts:
                pairs = [(a, b) for i, a in enumerate(small) for b in small[i:]]
                rnd.shuffle(pairs)
                for a, b in pairs[:400]:
                    t1 = v * a * b; e1 = -1 if t1 in T else 1
                    r2 = (r + e1 * contrib(v, t1)) % (v*v)
                    for o in level1_options(r2, excl=(t1,)):
                        if o[1]: opts.append((e1 + o[0], [(t1, e1)] + o[1]))
                    if len(opts) > 30: break
        else:
            # first fix the level-0 part with a subset of the terms v^2*a (a < v), each toggled
            cands = small if len(small) <= 12 else small[:12]
            for mask in range(1, 1 << len(cands)):
                ops = []; r2 = r; d = 0
                for i, a in enumerate(cands):
                    if mask >> i & 1:
                        t = v * v * a; eps = -1 if t in T else 1
                        ops.append((t, eps)); d += eps
                        r2 = (r2 + eps * pow(a, -1, v*v)) % (v*v)
                if r2 % v: continue
                for o in level1_options(r2, excl=tuple(t for t, _ in ops))[:20]:
                    opts.append((d + o[0], ops + o[1]))
        if not opts: return None
        def score(o):
            ds = sum(e / t for t, e in o[1])
            return o[0] - lam * ds + 0.6 * rnd.random()
        best = min(opts, key=score)
        for t, e in best[1]: tog(t, e); S += e / t
    # closing DP over R_p mod p^2 for p in F
    target = crt([(-R.get(p, 0)) % (p*p) for p in F])
    groups = []
    Fterms = [t for t in range(8, 7*7*5 + 1) if all(p in Fset for p in factorint(t)) and sum(factorint(t).values()) == 3 and len(factorint(t)) > 1]
    for t in Fterms:
        s = -1 if t in T else 1
        vec = [((s * contrib(p, t)) % (p*p)) if t % p == 0 else 0 for p in F]
        groups.append([(crt(vec), s, [(t, s)])])
    cof = [a*b for a, b in combinations_with_replacement(F, 2)]   # m with Omega(m) = 2 from F
    for v in primerange(F[-1] + 1, Vmax):
        opts = []
        for k in (1, 2, 3):
            for Q in combinations(cof, k):
                sg = [(-1 if v*m in T else 1) for m in Q]
                if sum(s_ * pow(m % v, -1, v) for s_, m in zip(sg, Q)) % v: continue
                vec = [0]*len(F)
                for s_, m in zip(sg, Q):
                    t = v * m
                    for i, p in enumerate(F):
                        if t % p == 0: vec[i] = (vec[i] + s_ * contrib(p, t)) % (p*p)
                opts.append((crt(vec), sum(sg), [(v*m, s_) for s_, m in zip(sg, Q)]))
        if opts: groups.append(opts)
    INF = 1e18; MU = 1.0 * X0
    dp = np.full(M, INF); dp[0] = 0.0; choice = []
    for g in groups:
        new = dp.copy(); ch = np.full(M, -1, dtype=np.int32)
        for oi, (e, c, ops) in enumerate(g):
            rolled = np.roll(dp, e)
            cand = np.where(rolled >= INF, INF, rolled + c - MU * sum(x / t for t, x in ops))
            better = cand < new; new[better] = cand[better]; ch[better] = oi
        choice.append(ch); dp = new
    if dp[target] >= INF / 2: return ('noclose', len(T), S)
    s = target; picks = []
    for gi in range(len(groups) - 1, -1, -1):
        oi = int(choice[gi][s])
        if oi >= 0:
            e, c, ops = groups[gi][oi]; picks.append(ops); s = (s - e) % M
    assert s == 0
    for ops in picks:
        for t, x in ops: tog(t, x)
    tot = sum(Fraction(1, t) for t in T)
    return (tot, len(T), T, S)
if __name__ == '__main__':
    X0, seed, lm = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    out = run(X0, seed, lm)
    if out is None: print('repair failed')
    elif out[0] == 'noclose': print('no close, terms %d, pre-sum %.4f' % (out[1], out[2]))
    else:
        tot, n, T, S0 = out
        print('sum', tot, 'terms', n, 'pre-sum %.4f' % S0, 'largest', max(T))
        if tot == 1 and len(sys.argv) > 4: open(sys.argv[4], 'w').write(' '.join(map(str, sorted(T))) + '\n')
