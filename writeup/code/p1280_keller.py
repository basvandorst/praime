"""Puzzle 1280: Keller cycles N -> (2N + P) mod 2P of Sierpinski, Riesel and Brier numbers.

Checks, with plain integers only:
  * every test number is covered by its covering set (one full period of n);
  * the Keller cycle has exactly M = lcm ord_p(2) members, all odd and below 2P,
    it closes at N mod 2P, and every member is covered by the same set;
  * Q1: the stairway N + jP (j = 2, 4, ...) never meets the cycle; the only cycle
    member in the class of N mod 2P is N mod 2P itself;
  * Q2: the mirror 2P - c of the cycle is the cycle of 2P - N and is covered with the
    opposite sign; Nmax - Nmin = 2P - Smin - Rmin; the gap formula through P/2;
  * all cycles of the covering sets {3,5,7,13,17,241} and {3,5,7,13,19,37,73};
  * every Brier number in Keller's catalogue (references/Catalogue.txt).
Run: python3 p1280_keller.py
"""
import os, re
from itertools import product
from math import gcd

def order2(p):
    # multiplicative order of 2 mod an odd prime p, by trial division of p - 1
    n, f, q = p - 1, [], 2
    while q * q <= n:
        if n % q == 0:
            f.append(q)
            while n % q == 0: n //= q
        q += 1
    if n > 1: f.append(n)
    t = p - 1
    for q in f:
        while t % q == 0 and pow(2, t // q, p) == 1: t //= q
    return t

def lcm(*a):
    r = 1
    for x in a: r = r * x // gcd(r, x)
    return r

def prod(a):
    r = 1
    for x in a: r *= x
    return r

def covered(k, cs, sign, M):
    return all(any((k * pow(2, n, p) + sign) % p == 0 for p in cs) for n in range(M))

def keller(N, P, M):
    x, c = N, []
    for _ in range(M):
        x = (2 * x + P) % (2 * P); c.append(x)
    return c

def check(N, sides, label, quiet=False):
    """sides: list of (covering set, sign); sign +1 Sierpinski, -1 Riesel."""
    CS = sorted(set(p for cs, s in sides for p in cs))
    P, M = prod(CS), lcm(*[order2(p) for p in CS])
    for cs, s in sides:
        Ms = lcm(*[order2(p) for p in cs])
        assert covered(N, cs, s, Ms), (label, cs, s)
        assert all(N % p for p in cs)                 # every prime of a minimal set is used
    C = keller(N, P, M)
    assert len(set(C)) == M and all(c % 2 == 1 and 0 < c < 2 * P for c in C)
    assert C[-1] == N % (2 * P)
    assert all(keller(N, P, t)[-1] != N % (2 * P) for t in range(1, M) if M % t == 0)
    for c in C:                                       # every member is covered
        for cs, s in sides:
            Ms = lcm(*[order2(p) for p in cs])
            assert covered(c, cs, s, Ms), (label, c)
    # Q1: stairway N + jP, j = 2, 4, ..., 40, against the cycle
    stair = [N + j * P for j in range(2, 41, 2)]
    assert all(t % (2 * P) == N % (2 * P) for t in stair)
    assert not set(stair) & set(C)
    assert [c for c in C if (c - N) % (2 * P) == 0] == [N % (2 * P)]
    # Q2: mirror
    mirror = [2 * P - c for c in C]
    assert set(mirror) == set(keller(2 * P - N % (2 * P), P, M))
    for c in mirror[:4]:
        for cs, s in sides:
            Ms = lcm(*[order2(p) for p in cs])
            assert covered(c, cs, -s, Ms)
    Smin, Nmax, Rmin = min(C), max(C), min(mirror)
    assert Nmax - Smin == 2 * P - Smin - Rmin
    # gap formula: orbit O = {2^i N mod P}; a = max of O below P/2, b = min of O above P/2
    O = {(N * pow(2, i, P)) % P for i in range(M)}
    a = max(o for o in O if 2 * o < P); b = min(o for o in O if 2 * o > P)
    assert Smin == 2 * b - P and Nmax == 2 * a + P and Smin + Rmin == 2 * (b - a)
    if not quiet:
        print(f"{label}: M = {M}, P = {P}, 2P = {2*P}")
        print(f"   cycle min {Smin}, max {Nmax}, range {Nmax - Smin} = 2P - ({Smin} + {Rmin}),"
              f" range/2P = {(Nmax - Smin)/(2*P):.6f}")
    return P, M, C, Smin, Nmax, Rmin

print("== Q1 and Q2 on single covering sets ==")
E6 = [3, 5, 7, 13, 17, 241]; S7 = [3, 5, 7, 13, 19, 37, 73]
check(78557, [(S7, 1)], "78557 (Sierpinski)")
check(271129, [(E6, 1)], "271129 (Sierpinski)")
check(509203, [(E6, -1)], "509203 (Riesel)")
check(19964725987718971, [([3, 5, 11, 13, 31, 41, 61, 151, 331, 1321], 1)], "19964725987718971 (puzzle 1279)")
check(11679013547, [([3, 5, 7, 11, 17, 31, 41, 61, 241, 1321], 1)], "11679013547 (A206001)")
check(10909924841, [([3, 5, 7, 13, 17, 19, 109, 193, 577, 1153], -1)], "10909924841 (A206430)")
B = 3316923598096294713661
SB = [3, 5, 13, 17, 97, 241, 673]; RB = [3, 7, 11, 19, 31, 37, 41, 73, 109, 151, 331, 1321]
P1, M1, *_ = check(B, [(SB, 1)], "Clavier 1, Sierpinski side")
P2, M2, *_ = check(B, [(RB, -1)], "Clavier 1, Riesel side")
P3, M3, *_ = check(B, [(SB, 1), (RB, -1)], "Clavier 1, Brier (both sides)")
assert (M1, P1) == (48, 52153970115) and (M2, P2) == (180, 108435121946131097271)
assert (M3, P3) == (720, 1885107369798300628981297352055)   # the page's numbers

print("\n== all Keller cycles of two covering sets ==")
def all_cycles(CS, sign):
    P = prod(CS); ms = [order2(p) for p in CS]; M = lcm(*ms)
    res = set()
    for a in product(*[range(m) for m in ms]):       # one class a_p (mod m_p) per prime
        cov = set()
        for m, ai in zip(ms, a): cov.update(range(ai, M, m))
        if len(cov) < M: continue
        N = 0                                         # CRT: N = -sign * 2^(-a_p) (mod p)
        for p, ai in zip(CS, a):
            r = (-sign * pow(2, -ai, p)) % p; q = P // p
            N = (N + r * q * pow(q, -1, p)) % P
        res.add(N if N % 2 else N + P)
    for N in res: assert covered(N, CS, sign, M)
    cycles, seen = [], set()
    for N in sorted(res):
        if N in seen: continue
        C = keller(N, P, M); seen.update(C); cycles.append(C)
    assert seen == res
    return P, M, res, cycles
for CS in (E6, S7):
    P, M, res, cyc = all_cycles(CS, 1)
    _, _, resR, cycR = all_cycles(CS, -1)
    assert {2 * P - x for x in res} == resR            # conversion R = 2P - S is a bijection
    print(f"{CS}: P = {P}, M = {M}, {len(res)} Sierpinski residues mod 2P in {len(cyc)} cycles")
    rows = []
    for C in cyc:
        _, _, _, Smin, Nmax, Rmin = check(min(C), [(CS, 1)], "", quiet=True)
        rows.append((Nmax - Smin, Smin, Nmax, Rmin))
    for r in sorted(rows, reverse=True):
        print(f"   Smin {r[1]:>9}  Nmax {r[2]:>9}  range {r[0]:>9} = 2P - {2*P - r[0]:>8}"
              f"  (Rmin of the mirror cycle {r[3]})  range/2P = {r[0]/(2*P):.5f}")
    print(f"   least Sierpinski residue {min(res)}, least Riesel residue {min(resR)}")

print("\n== Keller's catalogue of Brier numbers ==")
here = os.path.dirname(os.path.abspath(__file__))
lines = open(os.path.join(here, '..', '..', 'references', 'Catalogue.txt')).read().split('\n')
ents, name = [], None
for l in lines:
    s = l.strip()
    if s.startswith('k='):
        ents.append({'name': name, 'k': int(s[2:].split()[0]), 'cs': [], 'P': None})
    elif s.startswith('{') and ents:
        ents[-1]['cs'].append([int(x) for x in re.findall(r'\d+', s[:s.index('}')])])
    elif s.startswith('P=') and ents:
        ents[-1]['P'] = int(s[2:].split()[0])
    elif s and s[0].isalpha() and not s.startswith('Raw'):
        name = s
slips, stats, cmins = [], [], {}
for e in ents:
    if len(e['cs']) != 2: continue
    k, (A, Bc) = e['k'], e['cs']
    MA, MB = lcm(*map(order2, A)), lcm(*map(order2, Bc))
    if covered(k, A, 1, MA) and covered(k, Bc, -1, MB): S, R = A, Bc
    elif covered(k, Bc, 1, MB) and covered(k, A, -1, MA): S, R = Bc, A
    else:
        # Gallot 6: the listed Riesel set lacks 7 and 73
        assert e['name'] == 'Gallot 6' and covered(k, A, 1, MA) and not covered(k, Bc, -1, MB)
        R = sorted(Bc + [7, 73]); S = A
        assert covered(k, R, -1, lcm(*map(order2, R)))
        slips.append(f"{e['name']}: listed Riesel set misses 7 and 73; {{{','.join(map(str, R))}}} covers")
    P, M, C, Smin, Nmax, Rmin = check(k, [(S, 1), (R, -1)], e['name'], quiet=True)
    if P != e['P']:
        slips.append(f"{e['name']}: listed P is not the product of the listed primes (should be {P})")
    assert Smin == k                                    # the catalogue lists cycle minima
    stats.append((e['name'], M, P, Smin, Nmax, Rmin)); cmins[Smin] = e['name']
print(f"{len(stats)} Brier entries checked; each is the least member of its Keller cycle;"
      f" M from {min(r[1] for r in stats)} to {max(r[1] for r in stats)}")
basic = [r for r in stats if ' = ' not in r[0]]
assert len(basic) == 283
nk, nkP = len({r[3] for r in basic}), len({(r[3], r[2]) for r in basic})
print(f"283 basic entries: {nk} distinct k, {nkP} distinct (k, P); a k listed twice has two covering pairs")
for r in basic:                                       # no two basic numbers in one class
    C = set(keller(r[3], r[2], r[1])); Cm = {2 * r[2] - c for c in C}
    assert not [q for q in basic if q[3] != r[3] and q[2] == r[2] and (q[3] in C or q[3] in Cm)]
print("no two distinct basic k with the same P share a cycle or a mirrored cycle")
bymin = {r[0]: r for r in stats}
kof = {e['name']: e['k'] for e in ents}
for e in ents:                                        # the catalogue's duplicates
    if ' = ' not in e['name']: continue
    a, b = e['name'].split(' = '); o = bymin[b]; P, M = o[2], o[1]
    C = set(keller(o[3], P, M)); kd = e['k'] % (2 * P)
    rel = 'same number' if e['k'] == o[3] else 'in the cycle' if kd in C else \
          'in the mirrored cycle' if 2 * P - kd in C else None
    assert rel and e['P'] in (P, None)
    print(f"   duplicate {a}: {rel} of {b}")
print("catalogue slips:"); [print("  ", s) for s in slips]
pairs = [(r[0], cmins[r[5]]) for r in stats if r[5] in cmins]
print("entries whose mirror cycle is also listed:", pairs)
rr = [(r[4] - r[3]) / (2 * r[2]) for r in stats]
print(f"range/2P over the catalogue: min {min(rr):.6f}, max {max(rr):.7f}")
z = sorted(r[5] * (r[1] + 1) / (2 * r[2]) for r in stats)
print(f"Rmin*(M+1)/2P: mean {sum(z)/len(z):.3f}, median {z[len(z)//2]:.3f} (random model: 1 and ln 2 = 0.693)")
print("all checks passed")
