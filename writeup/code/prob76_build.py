# Problem 76: distinct odd semiprimes whose reciprocals sum to exactly 1, with as few terms as possible
# (writeup/prob-76.md).  View the semiprimes pq as edges of a graph on the odd primes.  The sum has no prime p
# in its denominator iff sum_{q ~ p} q^-1 = 0 (mod p) ("p is balanced"); if every vertex is balanced the sum is
# an integer, and a sum in (0, 2) is then exactly 1.
#   1. bulk: all odd squarefree semiprimes pq <= X0;
#   2. repair every vertex v > 13 from the largest down, by toggling edges to smaller primes or an "upward gadget"
#      (edge (v,u) with u > v plus edges (u,f), f in F = {3,5,7,11,13}, chosen so that u stays balanced),
#      scored by (change in edge count) - LAM*X0*(change in sum);
#   3. close F exactly: a dynamic programme over the 15015 residue states of F, using toggles of the edges
#      inside F and gadgets (v, Q): toggle the edges (v,q), q in Q subset of F, with sum s_q q^-1 = 0 (mod v).
# Usage: python3 prob76_build.py X0 SEED LAM [OUTFILE]      e.g.  python3 prob76_build.py 4100 8 0.7
# The run is deterministic for given (X0, SEED, LAM); it prints the number of terms or 'no close'.
import sys, random
from sympy import primerange, isprime
from fractions import Fraction
from itertools import combinations
from math import prod
F = [3, 5, 7, 11, 13]
M = prod(F)
def crt_vec(vec):
    s = 0
    for r, f in zip(vec, F):
        Mf = M // f
        s = (s + r * Mf * pow(Mf, -1, f)) % M
    return s
UMAX = 200000
Fset = set(F)
def run(X0, seed, lam_mult=1.0, Vfresh=200000):
    rnd = random.Random(seed)
    lam = lam_mult * X0
    P = list(primerange(3, X0 // 3 + 1))
    E = set(); nb = {p: set() for p in P}
    for i, p in enumerate(P):
        for q in P[i+1:]:
            if p*q > X0: break
            E.add((p, q)); nb[p].add(q); nb[q].add(p)
    def tog(v, u, s):
        e = (min(u, v), max(u, v))
        nb.setdefault(v, set()); nb.setdefault(u, set())
        if s == '+': E.add(e); nb[v].add(u); nb[u].add(v)
        else: E.discard(e); nb[v].discard(u); nb[u].discard(v)
    for v in sorted((p for p in P if p > F[-1]), reverse=True):
        if not nb[v]: continue
        r = sum(pow(u, -1, v) for u in nb[v]) % v
        if r == 0: continue
        small = [u for u in P if u < v]
        opts = []
        def single(r_, excl):
            out = []
            u0 = pow(r_, -1, v)
            if u0 not in excl and u0 in nb[v]: out.append((-1, [('-', u0)]))
            u1 = pow((-r_) % v, -1, v)
            if u1 not in excl and u1 >= 3 and u1 % 2 and u1 < v and isprime(u1) and u1 not in nb[v]: out.append((1, [('+', u1)]))
            return out
        if len(small) <= 14:
            inv = [pow(u, -1, v) for u in small]
            for m in range(1, 1 << len(small)):
                rr = r; ops = []; dn = 0
                for k, u in enumerate(small):
                    if m >> k & 1:
                        if u in nb[v]: rr -= inv[k]; ops.append(('-', u)); dn -= 1
                        else: rr += inv[k]; ops.append(('+', u)); dn += 1
                if rr % v == 0: opts.append((dn, ops))
        else:
            opts += single(r, ())
            for u in small:
                if u in nb[v]: r2 = (r - pow(u, -1, v)) % v; f_ = ('-', u); dn = -1
                else: r2 = (r + pow(u, -1, v)) % v; f_ = ('+', u); dn = 1
                if r2 == 0: opts.append((dn, [f_])); continue
                for d2, o2 in single(r2, (u,)): opts.append((dn + d2, [f_] + o2))
        # upward gadgets: edge (v,u) with u > v plus edges (u,f), f in F, keeping u balanced
        for sv in (1, -1):
            c = pow((-r * sv) % v, -1, v)           # need sv * u^-1 = -r (mod v)  =>  u = c (mod v)
            u = c
            while u <= v: u += v
            tries = 0
            while u < UMAX and tries < 4000:
                tries += 1
                if u % 2 and u not in Fset and isprime(u):
                    pres_vu = u in nb.get(v, ())
                    if (sv == 1 and not pres_vu) or (sv == -1 and pres_vu):
                        presu = nb.get(u, set())
                        vi = pow(v, -1, u)
                        for k in range(1, len(F) + 1):
                            for Q in combinations(F, k):
                                sg = [(-1 if f in presu else 1) for f in Q]
                                if (sv * vi + sum(x * pow(f, -1, u) for x, f in zip(sg, Q))) % u == 0:
                                    ops = [('+' if sv == 1 else '-', u)] + [('G', (f, u, x)) for f, x in zip(Q, sg)]
                                    opts.append((sv + sum(sg), ops))
                u += v
        if not opts:
            for i, u in enumerate(small):
                for w in small[i+1:]:
                    r3 = r; ops = []; dn = 0
                    for x in (u, w):
                        if x in nb[v]: r3 = (r3 - pow(x, -1, v)) % v; ops.append(('-', x)); dn -= 1
                        else: r3 = (r3 + pow(x, -1, v)) % v; ops.append(('+', x)); dn += 1
                    if r3 == 0: opts.append((dn, ops)); continue
                    for d3, o3 in single(r3, (u, w)): opts.append((dn + d3, ops + o3))
                if len(opts) > 40: break
        if not opts: return None
        def score(o):
            ds = 0.0
            for s_, u in o[1]:
                if s_ == 'G': f, uu, x = u; ds += x/(f*uu)
                else: ds += (1 if s_ == '+' else -1)/(v*u)
            return o[0] - lam*ds + 0.6*rnd.random()
        best = min(opts, key=score)
        for s_, u in best[1]:
            if s_ == 'G': f, uu, x = u; tog(f, uu, '+' if x > 0 else '-')
            else: tog(v, u, s_)
    # closing DP over residues at F
    res = [sum(pow(q, -1, f) for q in nb[f]) % f for f in F]
    target = crt_vec([(-x) % f for x, f in zip(res, F)])
    groups = []   # each group: list of (effect, cost, ops)
    for i, a in enumerate(F):
        for b in F[i+1:]:
            s = -1 if b in nb[a] else 1
            eff = [0]*len(F); eff[i] = (s * pow(b, -1, a)) % a; eff[F.index(b)] = (s * pow(a, -1, b)) % b
            groups.append([(crt_vec(eff), s, [(a, b, s)])])
    cand_v = [v for v in primerange(F[-1] + 1, Vfresh)]
    for v in cand_v:
        present = nb.get(v, set())
        if not present and v <= X0 // 3:
            pass
        opts = []
        for k in range(1, len(F) + 1):
            for Q in combinations(range(len(F)), k):
                sgn = [(-1 if F[i] in present else 1) for i in Q]
                if not present and any(x < 0 for x in sgn): continue
                if sum(sg * pow(F[i], -1, v) for sg, i in zip(sgn, Q)) % v: continue
                vinv = [pow(v, -1, F[i]) for i in Q]
                eff = [0]*len(F)
                for sg, i, vi in zip(sgn, Q, vinv): eff[i] = (sg * vi) % F[i]
                opts.append((crt_vec(eff), sum(sgn), [(F[i], v, sg) for sg, i in zip(sgn, Q)]))
        if opts: groups.append(opts)
    INF = 10**9
    dp = [INF]*M; dp[0] = 0
    choice = []
    for g in groups:
        new = dp[:]; ch = [-1]*M
        for oi, (e, c, ops) in enumerate(g):
            for s0 in range(M):
                if dp[s0] < INF:
                    s1 = (s0 + e) % M
                    if dp[s0] + c < new[s1]: new[s1] = dp[s0] + c; ch[s1] = oi
        choice.append(ch); dp = new
    S0 = sum(1.0/(p*q) for p, q in E)
    if dp[target] >= INF: return ('noclose', len(E), S0)
    # reconstruct
    s = target; picks = []
    for gi in range(len(groups) - 1, -1, -1):
        oi = choice[gi][s]
        if oi >= 0:
            e, c, ops = groups[gi][oi]; picks.append(ops); s = (s - e) % M
    assert s == 0
    for ops in picks:
        for a, b, sg in ops: tog(a, b, '+' if sg > 0 else '-')
    tot = sum(Fraction(1, p*q) for p, q in E)
    return (tot, len(E), set(E), len(groups))
if __name__ == '__main__':
    X0, seed, lm = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    out = run(X0, seed, lm)
    if out is None:
        print('repair failed')
    elif out[0] == 'noclose':
        print('no close (pre-closing sum %.3f)' % out[2])
    else:
        tot, n, E, ng = out
        print('sum', tot, 'terms', n, 'largest', max(p * q for p, q in E))
        if tot == 1 and len(sys.argv) > 4:
            open(sys.argv[4], 'w').write(' '.join(str(p * q) for p, q in sorted(E, key=lambda e: e[0] * e[1])) + '\n')
