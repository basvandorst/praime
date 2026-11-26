# review 162, task 1b: own root counts. a_K(p) = (number of degree-1 primes above p) - 1 for p <= 8|D|, p not dividing D.
# Method (independent of the authors' polrootsmod/idealprimedec): for a monic quartic g with p not dividing disc(g), the number of roots
# of g mod p is the trace of the Frobenius matrix on F_p[x]/(g) (each factor F_{p^d} contributes trace 1 if d = 1, else 0):
#   tr = 1 + [x^1] r + [x^2] r^2 + [x^3] r^3 mod g,  r = x^p mod (g, p), vectorised over all primes at once (int64, p < 2^19).
# g = polredabs polynomial (not the authors' polynomial); for p | index(g) a Tschirnhaus transform with p not dividing its index;
# p = 2, 3: direct evaluation. Output: per-field weighted sums in several binnings, and local types at l <= 47 (idealprimedec).
import sys, time, math, numpy as np
import cypari2
pari = cypari2.Pari(); pari.allocatemem(200_000_000)
R = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/'
r2 = int(sys.argv[1]); nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for q in range(2, int(n**0.5) + 1):
        if s[q]: s[q * q::q] = False
    return np.nonzero(s)[0].astype(np.int64)
PR = sieve(8 * 60000)

def mulmod(a, b, c, p):
    """a, b: lists of 4 arrays (deg < 4); c = [c0, c1, c2, c3] of monic g; returns a*b mod (g, p)."""
    d = [np.zeros_like(p) for _ in range(7)]
    for i in range(4):
        for j in range(4):
            d[i + j] = d[i + j] + a[i] * b[j]
    d = [x % p for x in d]
    for k in (6, 5, 4):
        for i in range(4):
            d[k - 4 + i] = (d[k - 4 + i] - d[k] * c[i]) % p
    return d[:4]

def timesx(a, c, p):
    top = a[3]; return [(-top * c[0]) % p, (a[0] - top * c[1]) % p, (a[1] - top * c[2]) % p, (a[2] - top * c[3]) % p]

def nroots_trace(coef, p):
    """coef: integer coefficients [c0, c1, c2, c3] of monic quartic; p: array of primes >= 5 with p not dividing disc."""
    c = [np.int64(x) % p for x in coef]
    r = [np.ones_like(p), np.zeros_like(p), np.zeros_like(p), np.zeros_like(p)]
    for b in range(int(p.max()).bit_length() - 1, -1, -1):
        r = mulmod(r, r, c, p); rx = timesx(r, c, p); bit = ((p >> b) & 1).astype(bool)
        r = [np.where(bit, rx[i], r[i]) for i in range(4)]
    r2_ = mulmod(r, r, c, p); r3_ = mulmod(r2_, r, c, p)
    return (1 + r[1] + r2_[2] + r3_[3]) % p

def coeffs(g):
    v = [int(x) for x in pari.Vec(g)]          # leading first
    assert v[0] == 1 and len(v) == 5
    return [v[4], v[3], v[2], v[1]]

def nroots_small(coef, q):
    return sum(1 for x in range(q) if (x**4 + coef[3] * x**3 + coef[2] * x**2 + coef[1] * x + coef[0]) % q == 0)

FE96 = np.linspace(0, 2, 193); FE192 = np.linspace(0, 2, 385)
ylo, yhi = (7 / 24)**3, (47 / 24)**3
YE = np.linspace(ylo, yhi, 21); SE = np.linspace(math.sqrt(ylo), math.sqrt(yhi), 21)
LT = [int(x) for x in sieve(47)]

def binsum(x, edges, w, wa):
    j = np.searchsorted(edges, x, side='left') - 1; n = len(edges) - 1; ok = (j >= 0) & (j < n)
    return np.bincount(j[ok], weights=wa[ok], minlength=n), np.bincount(j[ok], weights=w[ok], minlength=n)

rows = [l.split() for l in open(R + 'fields_r%d.txt' % r2)]
rng = np.random.default_rng(162 + r2); order = np.arange(len(rows))
if nmax < len(rows): order = np.sort(rng.choice(len(rows), nmax, replace=False))
out = {k: [] for k in ('n96', 'd96', 'n192', 'd192', 'ny', 'dy', 'ns', 'ds')}; Ds = []; polys = []; types = []; nfix = 0; t0 = time.time()
for i in order:
    D = int(rows[i][0]); g = pari(rows[i][1]); A = abs(D); cf = coeffs(g)
    p = PR[PR <= 8 * A]; p = p[A % p != 0]
    a = np.full(len(p), -99, np.int64)
    ind = math.isqrt(abs(int(pari.poldisc(g))) // A); assert ind * ind * A == abs(int(pari.poldisc(g)))
    big = p >= 5
    good = big & ((ind % p) != 0)
    a[good] = nroots_trace(cf, p[good]) - 1
    nf = pari.nfinit(g)
    for k in np.nonzero(~good)[0]:
        q = int(p[k]); h = g; tries = 0
        if q < 5:   # 2, 3 can be common index divisors: degree-1 primes from idealprimedec
            a[k] = sum(1 for P in pari.idealprimedec(nf, q) if int(P[3]) == 1) - 1; continue
        while True:
            dh = abs(int(pari.poldisc(h))); ih = math.isqrt(dh // A)
            if ih % q != 0 and pari.poldegree(h) == 4 and int(pari.pollead(h)) == 1: break
            zk = pari.nfbasis(g); cs = [int(v) for v in np.random.default_rng(tries).integers(-3, 4, 4)]
            h = pari('charpoly(Mod(%s, %s))' % ('+'.join('(%d)*(%s)' % (cc, zk[ii]) for ii, cc in enumerate(cs)), g)); tries += 1; assert tries < 50
            if int(pari.issquarefree(h)) == 0 or pari.poldegree(h) != 4: h = g; continue
        a[k] = nroots_small(coeffs(h), q) - 1; nfix += (q >= 5)
    assert (a >= -1).all() and (a <= 3).all()
    pf = p.astype(float); w = np.log(pf); wa = w * math.sqrt(A) * a; y = pf / A; u = np.cbrt(y)
    for key, x, e in (('96', u, FE96), ('192', u, FE192), ('y', y, YE), ('s', np.sqrt(y), SE)):
        nn, dd = binsum(x, e, w, wa); out['n' + key].append(nn); out['d' + key].append(dd)
    ty = []
    for l in LT:
        dec = pari.idealprimedec(nf, l); ef = sorted([(int(P[2]), int(P[3])) for P in dec], key=lambda t: (-t[0], -t[1]))
        ty.append(' '.join('%d^%d' % (f_, e_) if e_ > 1 else '%d' % f_ for e_, f_ in ef))
    types.append(ty); Ds.append(D); polys.append(rows[i][1])
np.savez_compressed(R + 'mine_r%d.npz' % r2, D=np.array(Ds), idx=order, types=np.array(types), LT=np.array(LT),
                    **{k: np.array(v, np.float64) for k, v in out.items()})
print('r2=%d: %d fields, %d index primes >= 5 fixed by another O_K generator, %.1f s' % (r2, len(Ds), nfix, time.time() - t0))
