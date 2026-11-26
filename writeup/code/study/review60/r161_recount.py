"""Review of round 161, task 5 (+ spike mechanism). Own field lists (PARI nflist), own root counts (x^p mod (f, p) and the
discriminant's Legendre symbol, vectorised over primes; brute force for p < 50; PARI idealprimedec only for index divisors).
Compares per-field bin sums with the stored ones; also bin sums without the primes represented by the principal form
(4p = t^2 + |D| n^2 for D < 0, 4p = |D| n^2 - t^2 for D > 0, n >= 1), and the a_K values at those primes."""
import os, math, random
import numpy as np, cypari2
H = os.path.dirname(os.path.abspath(__file__))
SC = "/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/cubic"
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
EDGES = np.round(np.linspace(0.0, 3.0, 61), 10)

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)
PR = sieve(130000)

def mulmod(u, v, A, B, C, p):
    u0, u1, u2 = u; v0, v1, v2 = v
    c0 = (u0 * v0) % p; c1 = (u0 * v1 + u1 * v0) % p; c2 = (u0 * v2 + u1 * v1 + u2 * v0) % p
    c3 = (u1 * v2 + u2 * v1) % p; c4 = (u2 * v2) % p
    # x^3 = A x^2 + B x + C ; x^4 = (A^2 + B) x^2 + (AB + C) x + AC   (A = -a, B = -b, C = -c mod p)
    A2B = (A * A + B) % p; ABC = (A * B + C) % p; AC = (A * C) % p
    r2 = (c2 + c3 * A + c4 * A2B) % p; r1 = (c1 + c3 * B + c4 * ABC) % p; r0 = (c0 + c3 * C + c4 * AC) % p
    return (r0, r1, r2)

def powmod_scalar(x, e, p):
    r = np.ones_like(p); x = x % p; e = e.copy()
    while (e > 0).any():
        m = (e & 1) == 1; r = np.where(m, (r * x) % p, r); x = (x * x) % p; e >>= 1
    return r

def aK(coef, D, P):
    a, b, c = coef
    disc = a * a * b * b - 4 * b ** 3 - 4 * a ** 3 * c - 27 * c * c + 18 * a * b * c
    assert disc % abs(D) == 0; ind2 = disc // D; ind = math.isqrt(ind2); assert ind * ind == ind2
    p = P; A = (-a) % p; B = (-b) % p; C = (-c) % p
    res = (np.ones_like(p), np.zeros_like(p), np.zeros_like(p)); base = (np.zeros_like(p), np.ones_like(p), np.zeros_like(p))
    e = p.copy()
    while (e > 0).any():
        m = (e & 1) == 1
        prod = mulmod(res, base, A, B, C, p); res = tuple(np.where(m, prod[i], res[i]) for i in range(3))
        base = mulmod(base, base, A, B, C, p); e >>= 1
    split = (res[0] == 0) & (res[1] == 1 % p) & (res[2] == 0)
    leg = powmod_scalar(np.int64(disc) % p, (p - 1) // 2, p)
    roots = np.where(split, 3, np.where(leg == p - 1, 1, 0))
    small = p < 50
    for i in np.nonzero(small)[0]:
        q = int(p[i]); roots[i] = sum(1 for x in range(q) if (x ** 3 + a * x * x + b * x + c) % q == 0)
    out = roots - 1
    special = [i for i in range(len(p)) if disc % int(p[i]) == 0 and abs(D) % int(p[i]) != 0]
    return out, special, disc

def field_rows(sgn):
    fl = pari('nflist("S3", [20000, 39999], %d)' % (1 if sgn < 0 else 0))
    rows = []
    for f in fl:
        D = int(pari.nfdisc(f)); rows.append((D, str(f).replace(' ', '')))
    return sorted(set(rows), key=lambda r: (abs(r[0]), r[1]))

def stored_rows(sgn):
    rows = [l.split() for l in open(SC + "/fields_hold_%+d.txt" % sgn) if not l.startswith('#')]
    return sorted(set((int(D), f) for D, f in rows), key=lambda r: (abs(r[0]), r[1]))

out = {}
for sgn in (-1, 1):
    mine = field_rows(sgn); st = stored_rows(sgn)
    z = np.load(SC + "/cubic_hold_%+d.npz" % sgn)
    print("sign %+d: nflist now %d fields, stored list %d, identical: %s; stored D order matches: %s" % (
        sgn, len(mine), len(st), mine == st, all(int(z['D'][i]) == st[i][0] for i in range(len(st)))))
    import sys
    random.seed(161 + sgn); idx = sorted(random.sample(range(len(st)), 300)) if 'all' not in sys.argv else list(range(len(st)))
    n_ = len(idx); NUM = np.zeros((n_, 60)); DEN = np.zeros((n_, 60)); NUMx = np.zeros((n_, 60)); DENx = np.zeros((n_, 60)); NUMp = np.zeros((3, n_, 60))
    aprin = []; maxd = 0.0; nspecial = 0
    for r, i in enumerate(idx):
        D, fs = st[i]; f = pari(fs); co = [int(x) for x in pari.Vec(f)]
        assert co[0] == 1 and len(co) == 4
        A = abs(D); P = PR[PR <= 3 * A]
        a, special, disc = aK(co[1:], D, P)
        for j in special:  # index divisors: degree-1 primes of K above p
            q = int(P[j]); nf = pari.nfinit(f); a[j] = sum(1 for pr in pari.idealprimedec(nf, q) if int(pr[3]) == 1) - 1; nspecial += 1
        keep = (A % P) != 0; p = P[keep]; a = a[keep]; y = p / A
        j = np.searchsorted(EDGES, y, side='left') - 1; ok = (j >= 0) & (j < 60)
        lp = np.log(p.astype(float))
        NUM[r] = np.bincount(j[ok], lp[ok] * math.sqrt(A) * a[ok], 60); DEN[r] = np.bincount(j[ok], lp[ok], 60)
        # principal-form primes
        prin = np.zeros(len(p), bool)
        for n in range(1, 4):
            s = 4 * p - A * n * n if sgn < 0 else A * n * n - 4 * p
            g = s >= 0; t = np.zeros_like(s); t[g] = np.round(np.sqrt(s[g].astype(float))).astype(np.int64)
            hit = g & (t * t == s); prin |= hit
            okn = ok & hit; NUMp[n - 1, r] = np.bincount(j[okn], lp[okn] * math.sqrt(A) * a[okn], 60)
        aprin.extend(a[prin].tolist())
        okx = ok & ~prin
        NUMx[r] = np.bincount(j[okx], lp[okx] * math.sqrt(A) * a[okx], 60); DENx[r] = np.bincount(j[okx], lp[okx], 60)
        maxd = max(maxd, float(np.max(np.abs(NUM[r] - z['num'][i]))), float(np.max(np.abs(DEN[r] - z['den'][i]))))
    rel = np.abs(NUM.sum(0) / DEN.sum(0) - z['num'][idx].sum(0) / z['den'][idx].sum(0)).max()
    av = np.array(aprin)
    print("   %d fields:" % n_ + "  max |own - stored| per field and bin: %.2e (num, den); pooled bin means max diff %.2e; index-divisor primes %d" % (maxd, rel, nspecial))
    print("   principal-form primes: %d, a_K values: %s" % (len(av), {int(v): int((av == v).sum()) for v in np.unique(av)}))
    np.savez(H + "/r161_%s_%+d.npz" % ("all" if "all" in sys.argv else "sample", sgn), idx=np.array(idx), NUM=NUM, DEN=DEN, NUMx=NUMx, DENx=DENx, NUMp=NUMp, D=np.array([st[i][0] for i in idx]))
