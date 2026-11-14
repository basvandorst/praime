"""Independent check of an LP-certified chain (mpmath and sympy; shares only the error budget delta with the chain).

1. The certified primes equal sympy's list.
2. The identity at random sample points of the chain's stages: zero side in mpmath from the 28-digit zero
   strings (mpmath's own Bessel functions), prime side exact; |difference| <= delta.
3. Proof replay: for recorded decisions, the window (unknowns, samples, A, y) is rebuilt in mpmath and every
   node of every refutation tree is re-verified with exact comparisons:
   forcings and violated rows (A >= 0), Farkas vectors (min over the box of (A^T w).b > w.y + |w|.delta),
   and the tree itself (each split has both children, every node descends from the root).
Usage: lp_check.py chain.json [max_decisions]
"""
import json, math, random, sys
import mpmath as mp
import numpy as np
from sympy import primerange, isprime, factorint, perfect_power
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from cert import sample   # the rigorous budget under test

mp.mp.dps = 30
import os
ZS = [l.strip() for l in open(os.environ.get('ZFILE', os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../data/riemann/zeros_hp2419.txt'))) if l.strip()]
Zm = [mp.mpf(s) for s in ZS]
Gld = np.array(ZS, dtype=np.longdouble)
EPS = float(os.environ.get('EPS', '2e-16'))   # the listed zeros are within EPS of the true ones (rounded lists: half a unit of the last decimal)

def Kmp(v, a, beta, m):
    t = 1 - (v / a) ** 2
    return t ** (mp.mpf(m) / 2) * mp.besseli(m, beta * mp.sqrt(t)) if t > 0 else mp.mpf(0)

def Khmp(w, a, beta, m):
    nu = m + mp.mpf(1) / 2; s = abs(a * w)
    if s < beta:
        q = mp.sqrt(beta ** 2 - s ** 2); z = mp.besseli(nu, q) / q ** nu if q > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    else:
        y = mp.sqrt(s ** 2 - beta ** 2); z = mp.besselj(nu, y) / y ** nu if y > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    return a * mp.sqrt(2 * mp.pi) * beta ** m * z

def Pmp(x, a, beta, m):
    nu = m + mp.mpf(1) / 2; r = a / (2 * x); q = mp.sqrt(beta ** 2 + r ** 2)
    return a * mp.sqrt(2 * mp.pi) * beta ** m * mp.besseli(nu, q) / q ** nu

_lam = {}
def lam(n):
    if n not in _lam:
        f = factorint(n); _lam[n] = mp.log(next(iter(f))) if len(f) == 1 else mp.mpf(0)
    return _lam[n]

def zero_side(x, a, beta, m, T):
    x = mp.mpf(x); L = mp.log(x); s = mp.mpf(0)
    for g in Zm:
        if g > T: break
        s += mp.cos(g * L) * Khmp(g / x, a, beta, m)
    triv = mp.quad(lambda v: Kmp(v, a, beta, m) * mp.sqrt(x * mp.e ** (v / x)) / (x * ((x * mp.e ** (v / x)) ** 3 - x * mp.e ** (v / x))), [-a, 0, a])
    return Pmp(x, a, beta, m) - 2 / mp.sqrt(x) * s - mp.sqrt(x) * triv

def support(x, a):
    x = mp.mpf(x); return int(mp.floor(x * mp.e ** (-a / x))), int(mp.ceil(x * mp.e ** (a / x)))

def prime_side(x, a, beta, m, skip=()):
    lo, hi = support(x, a); xm = mp.mpf(x)
    return sum(lam(n) * mp.sqrt(xm / n) * Kmp(xm * mp.log(mp.mpf(n) / xm), a, beta, m) for n in range(max(2, lo), hi + 1) if n not in skip)

out = []
def say(s):
    print(s, flush=True); out.append(s)

J = json.load(open(sys.argv[1])); maxdec = int(sys.argv[2]) if len(sys.argv) > 2 else 3
X, start = J['X'], J.get('start', 3)
wheel = J.get('wheel', [2])
# 1
hi_dec = J['result']['stop'] or X
true = [p for p in primerange(max(3, start), hi_dec) if all(p % q for q in wheel)]
got = [p for p in J['primes'] if p < hi_dec]
say(f'1. {len(got)} certified primes in [{max(3, start)}, {hi_dec}); equal to sympy: {got == true}; stop: {J["result"]["stop"]}; zeros up to {J["result"].get("maxT", 0):.2f}')

# 2
random.seed(7)
for st in J['stages'][-3:]:
    end, a, beta, m, cc, k = st
    x = float(random.randint(max(int(a * 1.5), end // 2), end - 1))
    T = cc * beta * x / a
    zs = zero_side(x, a, beta, m, T); ps = prime_side(x, a, beta, m)
    _, d, _ = sample(x, Gld, T, a, beta, m, EPS)
    k0 = mp.besseli(m, beta)
    say(f'2. x = {x:.0f} (a = {a}, beta = {beta}, m = {m}, T = {T:.1f}): |zero side - prime side| / K0 = {mp.nstr(abs(zs - ps) / k0, 3)}, delta / K0 = {d / float(k0):.3e}, inside: {abs(zs - ps) <= d}')

# 3
P = json.load(open(sys.argv[1].replace('.json', '_proofs.json')))
keys = sorted(P, key=lambda c: -sum(len(p['nodes']) for p in P[c]['proofs']))[:maxdec]   # the largest trees
for key in keys:
    rec = P[key]; U = rec['U']; xs = rec['xs']; a, beta, m, cc, k, h = rec['params']; c = U[0]
    # unknowns rebuilt independently: numbers >= c coprime to the wheel, not perfect powers, not yet decided
    V = []; n = c
    while len(V) < k + 1:
        if all(n % q for q in wheel) and not perfect_power(n): V.append(n)
        n += 1
    wnext = V[-1]
    okU = V[:-1] == U
    okX = all(mp.mpf(x) * mp.e ** (-a / mp.mpf(x)) > 1 and abs(x * mp.log(mp.mpf(c) / x)) < a and x * mp.log(mp.mpf(wnext) / x) >= a for x in xs)
    k0 = mp.besseli(m, beta)
    A = []; y = []; d = []
    for x in xs:
        T = cc * beta * x / a
        zs = zero_side(x, a, beta, m, T)
        known = prime_side(x, a, beta, m, skip=set(U))
        xm = mp.mpf(x)
        y.append((zs - known) / k0)
        d.append(mp.mpf(sample(x, Gld, T, a, beta, m, EPS)[1]) / k0)
    # column u: ln(u) sqrt(x/u) K(x ln(u/x)) / K0  (b_u = 1 means u is prime)
    A = [[mp.log(u) * mp.sqrt(mp.mpf(x) / u) * Kmp(mp.mpf(x) * mp.log(mp.mpf(u) / x), a, beta, m) / k0 for u in U] for x in xs]
    n_rows = len(xs); kk = len(U)
    lo = [y[i] - d[i] for i in range(n_rows)]; hi = [y[i] + d[i] for i in range(n_rows)]
    truth = [1 if isprime(u) else 0 for u in U]
    def rowsum(i, st, which):
        return sum(A[i][j] for j in range(kk) if (st[j] == 1 if which == 'base' else st[j] == -1))
    bad = []
    certified = []
    for pr in rec['proofs']:
        j, val, root = pr['bit'], pr['value'], pr['root']
        dec = lambda s: [{'x': -1, '0': 0, '1': 1}[ch] for ch in s]
        r0 = dec(root)
        if r0[j] != 1 - val or any(r0[i] != certified[i] for i in range(j)) or any(r0[i] != -1 for i in range(j + 1, kk)):
            bad.append(f'bit {j}: root state does not fix exactly the earlier bits and 1 - value')
        states = {}
        children_needed = []
        for node in pr['nodes']:
            st = dec(node['state']); states[node['state']] = True
            closed = False
            for ev in node['events']:
                if ev[0] == 'force':
                    _, var, v, row = ev
                    base = rowsum(row, st, 'base'); tot = rowsum(row, st, 'free')
                    ok = base + A[row][var] > hi[row] if v == 0 else base + tot - A[row][var] < lo[row]
                    if st[var] != -1 or not ok: bad.append(f'bit {j}: forcing of {U[var]} to {v} by row {row} fails')
                    st[var] = v
                else:
                    _, row, side = ev
                    base = rowsum(row, st, 'base'); tot = rowsum(row, st, 'free')
                    ok = base > hi[row] if side == 'hi' else base + tot < lo[row]
                    if not ok: bad.append(f'bit {j}: row {row} is not violated')
                    closed = True
            if closed:
                continue
            if 'farkas' in node:
                w = [mp.mpf(v) for v in node['farkas']]
                free = [q for q in range(kk) if st[q] == -1]
                g = [sum(w[i] * A[i][q] for i in range(n_rows) if w[i] != 0) for q in free]
                minbox = sum(min(gq, 0) for gq in g)
                rhs = sum(w[i] * (y[i] - rowsum(i, st, 'base')) + abs(w[i]) * d[i] for i in range(n_rows) if w[i] != 0)
                if not minbox > rhs: bad.append(f'bit {j}: Farkas vector fails ({mp.nstr(minbox - rhs, 3)})')
            elif 'branch' in node:
                q = node['branch']
                for v in (0, 1):
                    s2 = list(st); s2[q] = v
                    children_needed.append(''.join('x01'[t + 1] for t in s2))
            else:
                bad.append(f'bit {j}: open node')
        missing = [s for s in children_needed if s not in states]
        if missing: bad.append(f'bit {j}: {len(missing)} children missing')
        if len(states) != 1 + len(children_needed): bad.append(f'bit {j}: node count {len(states)} != 1 + {len(children_needed)}')
        certified.append(val)
    agree = certified == truth[:len(certified)]
    say(f'3. window at {c} (a = {a}, beta = {beta}, m = {m}, {n_rows} samples, {kk} unknowns): unknowns rebuilt equal: {okU}; samples valid: {okX}; '
        f'{len(certified)} decisions replayed, {sum(len(p["nodes"]) for p in rec["proofs"])} nodes, failures: {bad[:3] if bad else 0}; decisions equal sympy: {agree}')
open(sys.argv[1].replace('.json', '_check.txt'), 'w').write('\n'.join(out) + '\n')
