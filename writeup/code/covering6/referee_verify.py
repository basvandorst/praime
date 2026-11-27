# Referee's own verifier of the Lemma 2 certificate (written from proof.txt only).
# usage: python3 -I myverify.py CERT [extra_moduli_comma_sep] [inverse]
import sys, json
from collections import Counter

L = 180
DIVS = [d for d in range(1, L + 1) if L % d == 0]
ODD = [1, 3, 5, 9, 15, 45]
T = [4 * d for d in ODD]
B = [d for d in DIVS if d >= 6]
extra = []
if len(sys.argv) > 2 and sys.argv[2] not in ('-', ''):
    extra = [int(x) for x in sys.argv[2].split(',')]
B = sorted(set(B) | set(extra))
USE_INV = len(sys.argv) > 3 and sys.argv[3] == 'inverse'
assert B[:13] == [6, 9, 10, 12, 15, 18, 20, 30, 36, 45, 60, 90, 180] or extra, B
assert T == [4, 12, 20, 36, 60, 180]

def isprime(q):
    return q >= 2 and all(q % p for p in range(2, q))

FACT = {4: 2, 9: 3, 5: 5}   # factor size -> prime
def crt(x4, x9, x5):
    for b in range(L):
        pass
CRT = {}
for b in range(L):
    CRT[(b % 4, b % 9, b % 5)] = b
assert len(CRT) == L

def make_phi(perm):
    n = len(perm)
    assert n in (4, 9, 5), n
    assert sorted(perm) == list(range(n)), perm
    pi = list(perm)
    if USE_INV:
        inv = [0] * n
        for x, y in enumerate(pi):
            inv[y] = x
        pi = inv
    phi = []
    for b in range(L):
        c = [b % 4, b % 9, b % 5]
        idx = {4: 0, 9: 1, 5: 2}[n]
        c[idx] = pi[c[idx]]
        phi.append(CRT[tuple(c)])
    assert sorted(phi) == list(range(L))
    return phi

def classes_preserved(phi):
    # phi maps every residue class mod every divisor of 180 onto a residue class of the same modulus
    for m in DIVS:
        for r in range(m):
            img = {phi[b] for b in range(L) if b % m == r}
            rs = {x % m for x in img}
            if len(rs) != 1:
                return False
    return True

def rset(g, rho):
    return frozenset(b for b in range(L) if b % g == rho % g)

stats = Counter()
phis = set()
maxw = 0
maxdepth = 0
fail = []
phi_ok_cache = {}

def leaf(state, W, path):
    global maxw
    for k, v in W.items():
        assert isinstance(v, int) and v >= 0, (k, v)
        b = int(k)
        assert 0 <= b < L
        maxw = max(maxw, v)
    w = [0] * L
    for k, v in W.items():
        w[int(k)] = v
    fixed = [m for m in B if state[m][0] == m]
    U = [b for b in range(L) if all(b % m != state[m][1] % m for m in fixed)]
    WU = sum(w[b] for b in U)
    S = 0
    for m in B:
        g, rho = state[m]
        if g == m:
            continue
        best = 0
        for r in range(m):
            if r % g != rho % g:
                continue
            s = sum(w[b] for b in U if b % m == r)
            best = max(best, s)
        S += best
    N = WU - S
    D = 0
    for t in T:
        D += max(sum(w[b] for b in U if b % t == r) for r in range(t))
    stats['leaf'] += 1
    if N == D:
        stats['tight'] += 1
    if not (D > 0 and N >= D):
        fail.append(('leaf', path, N, D))

def walk(node, state, path, depth):
    global maxdepth
    maxdepth = max(maxdepth, depth)
    if 'w' in node:
        assert set(node.keys()) == {'w'}, node.keys()
        leaf(state, node['w'], path)
        return
    assert 'sym' not in node
    assert set(node.keys()) == {'m', 'q', 'ch'}, node.keys()
    m, q, ch = node['m'], node['q'], node['ch']
    assert m in B and isprime(q), (m, q)
    g, rho = state[m]
    assert m % (g * q) == 0, (m, g, q)
    assert len(ch) == q, (len(ch), q)
    stats['int'] += 1
    kids = []
    for j in range(q):
        st = dict(state)
        st[m] = (g * q, (rho + g * j) % (g * q))
        kids.append(st)
    for j, c in enumerate(ch):
        if 'sym' in c:
            assert set(c.keys()) == {'sym'}
            j2, perm = c['sym']
            assert isinstance(j2, int) and 0 <= j2 < q and j2 != j
            assert 'sym' not in ch[j2], 'target is symmetric'
            phi = make_phi(perm)
            key = tuple(phi)
            if key not in phi_ok_cache:
                phi_ok_cache[key] = classes_preserved(phi)
            phis.add(key)
            ok = phi_ok_cache[key]
            # restriction sets: unbranched m must be invariant; branched m: child j set -> child j2 set
            for mm in B:
                Sj = rset(*kids[j][mm])
                Sj2 = rset(*kids[j2][mm])
                img = frozenset(phi[b] for b in Sj)
                if img != Sj2:
                    ok = False
                    break
                if mm != m:
                    assert Sj == Sj2
            stats['sym'] += 1
            if not ok:
                fail.append(('sym', path + [j], j2))
        else:
            walk(c, kids[j], path + [j], depth + 1)

cert = json.load(open(sys.argv[1]))
assert cert.get('n') == 2
state0 = {m: (1, 0) for m in B}
walk(cert['tree'], state0, [], 0)
print(dict(stats), 'distinct phi', len(phis), 'maxdepth', maxdepth, 'maxw', maxw)
print('B =', B, 'inverse perms' if USE_INV else '')
if fail:
    print('FAILURES', len(fail), 'first', fail[:3])
    print('REJECT')
else:
    print('ACCEPT')
