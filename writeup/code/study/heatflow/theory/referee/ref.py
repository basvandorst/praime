"""Referee checks for heatflow/theory/proof.md. Own code; nothing imported from the repo.
Sections: A Theorem 1 (exact, sympy), B circle property, C second moment (exact Weyl grid + MC),
D scaling limit K_tau(u) (own quadrature vs exact finite N)."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import sys, math, time, itertools
from fractions import Fraction
from functools import lru_cache
from collections import Counter
import numpy as np
import sympy as sp
import mpmath as mp

T0 = time.time()
NCHK = [0, 0]


def check(name, ok, info=''):
    NCHK[0] += 1
    if not ok:
        NCHK[1] += 1
    print(('PASS ' if ok else 'FAIL ') + name + ((' | ' + info) if info else ''), flush=True)


def partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def zee(lam):
    z = 1
    for i, m in Counter(lam).items():
        z *= i ** m * math.factorial(m)
    return z


def mfact(lam):
    z = 1
    for i, m in Counter(lam).items():
        z *= math.factorial(m)
    return z


# ---------------------------------------------------------------- graphs
def connected(r, edges, mask):
    adj = [0] * r
    for b, (s, s2) in enumerate(edges):
        if mask >> b & 1:
            adj[s] |= 1 << s2
            adj[s2] |= 1 << s
    seen, front = 1, 1
    while front:
        nf = 0
        for v in range(r):
            if front >> v & 1:
                nf |= adj[v]
        nf &= ~seen
        seen |= nf
        front = nf
    return seen == (1 << r) - 1


@lru_cache(None)
def conn_graphs(r):
    edges = list(itertools.combinations(range(r), 2))
    out = []
    for mask in range(1 << len(edges)):
        if connected(r, edges, mask):
            out.append(tuple(edges[b] for b in range(len(edges)) if mask >> b & 1))
    return tuple(out)


# polynomials in Q as dict exponent -> Fraction
def pmul(A, B):
    C = {}
    for a, x in A.items():
        for b, y in B.items():
            C[a + b] = C.get(a + b, 0) + x * y
    return {e: c for e, c in C.items() if c != 0}


def padd(A, B, s=1):
    C = dict(A)
    for b, y in B.items():
        C[b] = C.get(b, 0) + s * y
    return {e: c for e, c in C.items() if c != 0}


def phi_enum(lam):
    """literal enumeration of connected graphs, edge weight Q^{2 i i'} - 1"""
    r = len(lam)
    cnt = Counter()
    for G in conn_graphs(r):
        cnt[tuple(sorted(2 * lam[s] * lam[s2] for s, s2 in G))] += 1
    tot = {}
    for key, c in cnt.items():
        P = {0: Fraction(1)}
        for ex in key:
            P = pmul(P, {ex: Fraction(1), 0: Fraction(-1)})
        tot = padd(tot, {e: c * v for e, v in P.items()})
    return tot


def phi_rec(lam):
    """connected part by the rooted set recursion; F(S) = Q^{(sum_S i)^2 - sum_S i^2}"""
    r = len(lam)
    full = (1 << r) - 1

    def F(S):
        a = [lam[s] for s in range(r) if S >> s & 1]
        return {sum(a) ** 2 - sum(x * x for x in a): Fraction(1)}
    Phi = {}
    for S in range(1, full + 1):
        low = S & -S
        val = F(S)
        T = (S - 1) & S
        while T:
            if T & low and T != S:
                val = padd(val, pmul(Phi[T], F(S ^ T)), -1)
            T = (T - 1) & S
        Phi[S] = val
    return Phi[full]


def c_graph(lam, phi):
    r, k = len(lam), sum(lam)
    pre = Fraction((-1) ** (r - 1) * k, zee(lam))
    sh = sum(i * i for i in lam)
    return {e + sh: pre * v for e, v in phi.items()}


# ============================================================ A. Theorem 1
print('=== A. Theorem 1: exact symbolic check ===')
KMAX = 7
Qs = sp.Symbol('Q')
Ps = sp.symbols('p1:%d' % (KMAX + 1))
gens = list(Ps) + [Qs]
one = sp.Poly(1, *gens)
e = [one]
for j in range(1, KMAX + 1):   # Newton: j e_j = sum (-1)^{i-1} e_{j-i} p_i
    acc = sp.Poly(0, *gens)
    for i in range(1, j + 1):
        acc = acc + e[j - i] * sp.Poly((-1) ** (i - 1) * Ps[i - 1], *gens)
    e.append(acc * sp.Rational(1, j))
g = [e[j] * sp.Poly(Qs ** (j * j), *gens) for j in range(KMAX + 1)]
lg = [None] * (KMAX + 1)
for k in range(1, KMAX + 1):
    acc = sp.Poly(0, *gens)
    for j in range(1, k):
        acc = acc + lg[j] * g[k - j] * j
    lg[k] = g[k] - acc * sp.Rational(1, k)
qpoly = {k: lg[k] * ((-1) ** (k - 1) * k) for k in range(1, KMAX + 1)}

CLAM = {}   # exact c_lambda as Q-polys (log route)
for k in range(1, KMAX + 1):
    d = qpoly[k].as_dict()
    for lam in partitions(k):
        m = Counter(lam)
        mon = tuple(m.get(i, 0) for i in range(1, KMAX + 1))
        CLAM[lam] = {ex[-1]: Fraction(int(c.p), int(c.q)) for ex, c in d.items() if ex[:-1] == mon}
    # no other monomials
    allmons = set(ex[:-1] for ex in d)
    ok = all(sum((i + 1) * mm for i, mm in enumerate(mon)) == k for mon in allmons)
    check('q_%d is homogeneous of weight %d in the p-basis' % (k, k), ok)

counts = [len(conn_graphs(r)) for r in range(1, 7)]
check('connected graph counts r=1..6', counts == [1, 1, 4, 38, 728, 26704], str(counts))

nok = nlam = nbad_sym = nbad_sign = 0
for k in range(1, 7):
    for lam in partitions(k):
        nlam += 1
        cg = c_graph(lam, phi_enum(lam))
        nok += (cg == CLAM[lam])
        # sensitivity: drop the symmetry factor, or use sign (-1)^{k-1}
        wrong1 = {e_: v * mfact(lam) for e_, v in cg.items()}
        wrong2 = {e_: v * (-1) ** (k - len(lam)) for e_, v in cg.items()}
        nbad_sym += (wrong1 != CLAM[lam])
        nbad_sign += (wrong2 != CLAM[lam])
check('Theorem 1 (graph enumeration) == log route, all partitions k<=6, exact in Q', nok == nlam,
      '%d/%d; dropping 1/prod m_i! breaks %d, sign (-1)^(k-1) breaks %d' % (nok, nlam, nbad_sym, nbad_sign))
nok7 = sum(c_graph(lam, phi_rec(lam)) == CLAM[lam] for lam in partitions(7))
check('Theorem 1 for k=7 (rooted-set recursion for Phi) == log route', nok7 == 15, '%d/15' % nok7)
nrec = sum(phi_rec(lam) == phi_enum(lam) for k in range(1, 7) for lam in partitions(k))
check('rooted-set recursion == enumeration for Phi, k<=6', nrec == 29)
check('k=2: q_2 = Q^4 p_2 - (Q^4-Q^2) p_1^2', CLAM[(2,)] == {4: 1} and CLAM[(1, 1)] == {4: -1, 2: 1})

# numeric: k > N, complex z (not unimodular), complex t; p_k(t) from the ROOTS of E_t
mp.mp.dps = 40


def evalQ(P, Qv):
    return mp.fsum(mp.mpf(c.numerator) / c.denominator * Qv ** e_ for e_, c in P.items())


rng = np.random.default_rng(12345)
worst = 0
for N in (2, 3):
    for trial in range(2):
        z = [mp.mpc(float(a), float(b)) for a, b in rng.normal(size=(N, 2))]
        t = mp.mpc(0.07, -0.11) if trial == 0 else mp.mpc(-0.13, 0.05)
        el = [mp.mpc(1)] + [mp.fsum(mp.fprod(c) for c in itertools.combinations(z, j)) for j in range(1, N + 1)]
        coef = [mp.exp(-t * j * (N - j)) * el[j] for j in range(N + 1)]
        roots = mp.polyroots(coef[::-1], maxsteps=200, extraprec=200)
        zeta = [-1 / x for x in roots]
        Qv = mp.exp(t)
        for k in range(1, 7):
            lhs = mp.fsum(zz ** k for zz in zeta)
            rhs = mp.exp(-t * N * k) * mp.fsum(evalQ(CLAM[lam], Qv) * mp.fprod(mp.fsum(zz ** i for zz in z) for i in lam)
                                             for lam in partitions(k))
            worst = max(worst, float(abs(lhs - rhs) / abs(lhs)))
check('Theorem 1 numerically vs roots of E_t, N=2,3, k=1..6 (k>N incl.), complex z and t', worst < 1e-25,
      'max rel err %.1e' % worst)

# ============================================================ B. circle property
print('=== B. roots of E_t ===')


def haar(n, N, rng):
    Z = (rng.normal(size=(n, N, N)) + 1j * rng.normal(size=(n, N, N))) / np.sqrt(2)
    Qm, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2)
    return Qm * (d / np.abs(d))[:, None, :]


def elem(lamb):
    """elementary symmetric e_0..e_N of rows of lamb (n,N) -> (n,N+1)"""
    n, N = lamb.shape
    E = np.zeros((n, N + 1), complex)
    E[:, 0] = 1
    for i in range(N):
        E[:, 1:] = E[:, 1:] + lamb[:, i:i + 1] * E[:, :-1]
    return E


rng = np.random.default_rng(777)
U = haar(400, 7, rng)
ev = np.linalg.eigvals(U)
E7 = elem(ev)
N = 7
for t in (0.002, 0.1, 1.0):
    w = np.exp(-t * np.arange(N + 1) * (N - np.arange(N + 1)))
    dev = max(np.max(np.abs(np.abs(np.roots((E7[m] * w)[::-1])) - 1)) for m in range(400))
    check('t=%g: all roots on |x|=1 (400 Haar U, N=7)' % t, dev < 1e-9, 'max dev %.1e' % dev)
t = -0.05
w = np.exp(-t * np.arange(N + 1) * (N - np.arange(N + 1)))
off = sum(np.max(np.abs(np.abs(np.roots((E7[m] * w)[::-1])) - 1)) > 1e-6 for m in range(400))
check('t=-0.05: some Haar U have roots off the circle', off > 0, '%d/400 off' % off)
# gradient flow of log|Vandermonde|^2 by finite differences, one U, N=6, t0=0.25
U6 = haar(1, 6, np.random.default_rng(5))[0]
e6 = elem(np.linalg.eigvals(U6)[None, :])[0]
N = 6


def angs(t):
    w = np.exp(-t * np.arange(N + 1) * (N - np.arange(N + 1)))
    return np.sort(np.angle(np.roots((e6 * w)[::-1])))


t0, h = 0.25, 1e-5
th = angs(t0)
d_fd = (angs(t0 + h) - angs(t0 - h)) / (2 * h)
d_th = np.array([sum(1 / np.tan((th[i] - th[j]) / 2) for j in range(N) if j != i) for i in range(N)])
check('zero dynamics dtheta/dt = sum cot((theta_i-theta_j)/2)', np.max(np.abs(d_fd - d_th)) < 1e-6 * np.max(np.abs(d_th)),
      'max abs diff %.1e' % np.max(np.abs(d_fd - d_th)))

# ============================================================ C. second moment
print('=== C. Theorem 3 (second moment) ===')


@lru_cache(None)
def chi(nu, lam):
    if not lam:
        return 1 if not nu else 0
    h, rest = lam[0], lam[1:]
    L = len(nu)
    beta = [nu[i] + (L - 1 - i) for i in range(L)]
    bs = set(beta)
    tot = 0
    for b in beta:
        nb = b - h
        if nb >= 0 and nb not in bs:
            sign = (-1) ** sum(1 for c in beta if nb < c < b)
            newb = sorted([c for c in beta if c != b] + [nb], reverse=True)
            Ln = len(newb)
            newnu = tuple(x for x in (newb[i] - (Ln - 1 - i) for i in range(Ln)) if x > 0)
            tot += sign * chi(newnu, rest)
    return tot


okc = True
for k in range(1, 8):
    P = list(partitions(k))
    for a in P:
        for b in P:
            s = sum(chi(nu, a) * chi(nu, b) for nu in P)
            okc &= (s == (zee(a) if a == b else 0))
check('Murnaghan-Nakayama characters: column orthogonality S_1..S_7', okc)


def moments_formula(N, t, kmax):
    """returns dict k -> (e^{-2ReT N k} * [truncated Schur sum], e^{-2ReT N k} H_k)"""
    Qv = mp.exp(t)
    out = {}
    for k in range(1, kmax + 1):
        P = list(partitions(k))
        c = {lam: evalQ(CLAM[lam], Qv) for lam in P}
        full = mp.fsum(abs(c[l]) ** 2 * zee(l) for l in P)
        trunc = mp.fsum(abs(mp.fsum(c[l] * chi(nu, l) for l in P)) ** 2 for nu in P if len(nu) <= N)
        f = mp.exp(-2 * mp.re(t) * N * k)
        out[k] = (f * trunc, f * full)
    return out


def weyl_grid(N, tlist, kmax, M):
    """exact Haar average of |p_k(t)|^2 by Weyl's formula on the M-th roots of unity grid"""
    roots = np.exp(2j * np.pi * np.arange(M) / M)
    acc = {(ti, k): 0.0 for ti in range(len(tlist)) for k in range(1, kmax + 1)}
    norm = 0.0
    rest = (np.array(list(itertools.product(range(M), repeat=N - 1)), dtype=int) if N > 1
            else np.zeros((1, 0), dtype=int))
    for i0 in range(M):
        idx = np.column_stack([np.full(len(rest), i0), rest])
        z = roots[idx]
        V = np.ones(len(z))
        for a in range(N):
            for b in range(a + 1, N):
                V *= np.abs(z[:, a] - z[:, b]) ** 2
        E = elem(z)
        norm += V.sum()
        for ti, t in enumerate(tlist):
            jj = np.arange(N + 1)
            a = E * np.exp(-t * jj * (N - jj))
            a = np.concatenate([a, np.zeros((len(z), kmax + 1))], axis=1)
            l = [None] * (kmax + 1)
            for k in range(1, kmax + 1):
                s = a[:, k].copy()
                for j in range(1, k):
                    s -= j * l[j] * a[:, k - j] / k
                l[k] = s
                pk = (-1) ** (k - 1) * k * l[k]
                acc[(ti, k)] += np.sum(np.abs(pk) ** 2 * V)
    tot = M ** N * math.factorial(N)
    return {key: v / tot for key, v in acc.items()}, norm / tot


N5, t5 = 5, 0.02
tc = 0.02 + 0.03j
tg = time.time()
grid, nrm = weyl_grid(N5, [t5, tc], 7, 23)
print('  Weyl grid N=5, M=23: %.1f s, normalisation %.15f' % (time.time() - tg, nrm))
check('Weyl grid normalisation = 1', abs(nrm - 1) < 1e-12, '%.2e' % (nrm - 1))
form = moments_formula(N5, mp.mpf(t5), 7)
formc = moments_formula(N5, mp.mpc(tc.real, tc.imag), 5)

# Monte Carlo, 100k Haar U(5)
rng = np.random.default_rng(20261127)
nmc = 100000
vals = {k: [] for k in range(1, 8)}
for chunk in range(4):
    Ub = haar(nmc // 4, N5, rng)
    Eb = elem(np.linalg.eigvals(Ub))
    jj = np.arange(N5 + 1)
    a = np.concatenate([Eb * np.exp(-t5 * jj * (N5 - jj)), np.zeros((len(Eb), 8))], axis=1)
    l = [None] * 8
    for k in range(1, 8):
        s = a[:, k].copy()
        for j in range(1, k):
            s -= j * l[j] * a[:, k - j] / k
        l[k] = s
        vals[k].append(np.abs((-1) ** (k - 1) * k * l[k]) ** 2)
print('  N=5 t=0.02:  k | MC +- SE | Weyl grid | formula (trunc for k>N) | full Hall norm | z(MC) | r=1 term')
for k in range(1, 8):
    v = np.concatenate(vals[k])
    mc, se = v.mean(), v.std() / np.sqrt(len(v))
    tr, fu = float(form[k][0]), float(form[k][1])
    gr = grid[(0, k)]
    r1 = k * math.exp(-2 * t5 * k * (N5 - k)) if k <= N5 else float('nan')
    print('   %d | %.6f +- %.6f | %.10f | %.10f | %.10f | %+.2f | %.6f' % (k, mc, se, gr, tr, fu, (mc - tr) / se, r1))
    check('k=%d: Weyl grid == formula (%s)' % (k, 'H_k' if k <= N5 else 'truncated'), abs(gr - tr) < 1e-10 * tr,
          'rel %.1e' % (abs(gr - tr) / tr))
    check('k=%d: MC within 4 SE of formula' % k, abs(mc - tr) < 4 * se, 'z=%.2f' % ((mc - tr) / se))
    if k <= N5:
        check('k=%d<=N: truncated == full Hall norm' % k, abs(form[k][0] - form[k][1]) < mp.mpf('1e-30') * form[k][1])
    else:
        check('k=%d>N: full Hall norm differs from Haar value' % k, abs(fu - gr) > 1e-3 * gr,
              'full %.6f vs %.6f; full rejected by MC at z=%.0f' % (fu, gr, (mc - fu) / se))
for k in (3, 5):
    gr, fo = grid[(1, k)], float(formc[k][1])
    check('complex t=0.02+0.03i, k=%d: e^{-2Re t Nk} sum |c|^2 z == Weyl grid' % k, abs(gr - fo) < 1e-10 * fo,
          '%.10f vs %.10f' % (gr, fo))

# author's exact example N=4, k=6, Q=11/10 (values of E|q_6|^2 and H_6)
Qr = Fraction(11, 10)
P6 = list(partitions(6))
c6 = {lam: sum(c * Qr ** e_ for e_, c in CLAM[lam].items()) for lam in P6}
H6 = sum(c6[l] ** 2 * zee(l) for l in P6)
T6 = sum(sum(c6[l] * chi(nu, l) for l in P6) ** 2 for nu in P6 if len(nu) <= 4)
check('author example N=4,k=6,Q=11/10: truncated 428.570651, full 14529.767197',
      abs(float(T6) - 428.570651) < 1e-6 and abs(float(H6) - 14529.767197) < 1e-6, '%.6f, %.6f' % (float(T6), float(H6)))

# small-t expansion H_k = k e^{2tk^2} + t^2 k^2 (k^3-k)/3 + O(t^3), k=5
mp.mp.dps = 50
P5 = list(partitions(5))
r_ = []
for tt in (mp.mpf('2e-5'), mp.mpf('1e-5')):
    Hh = mp.fsum(evalQ(CLAM[l], mp.exp(tt)) ** 2 * zee(l) for l in P5)
    r_.append((Hh - 5 * mp.exp(2 * tt * 25)) / tt ** 2)
check('small t: (H_5 - 5e^{50t})/t^2 -> 25*120/3 = 1000', abs(2 * r_[1] - r_[0] - 1000) < 1e-3,
      '%.6f, %.6f, Richardson %.6f' % (r_[0], r_[1], 2 * r_[1] - r_[0]))
# bound H_k <= k^2 e^{2tk^2} (t>=0)
for (tt, k) in ((0.3, 6), (0.5, 5), (1.0, 4), (0.05, 7)):
    Pk = list(partitions(k))
    Hh = mp.fsum(evalQ(CLAM[l], mp.exp(tt)) ** 2 * zee(l) for l in Pk)
    ratio = float(Hh / (k * k * mp.exp(2 * tt * k * k)))
    check('H_%d(%g) <= k^2 e^{2tk^2}' % (k, tt), ratio <= 1, 'ratio %.4f' % ratio)

# "strictly larger for k > N": counterexamples at t<0
print('  -- k>N: is the full Hall norm always strictly larger? --')
mp.mp.dps = 40
xs = sorted(float(x) for x in np.roots([1, 0, -1, 1 / 3]) if abs(x.imag) < 1e-12 and x.real > 0)
for N, k, Q2 in [(1, 2, 0.5)] + [(2, 3, x) for x in xs]:
    tstar = mp.log(mp.mpf(Q2)) / 2 if N == 1 else mp.log(mp.findroot(lambda x: x ** 3 - x + mp.mpf(1) / 3, Q2)) / 2
    res = []
    for tt in (tstar, tstar + mp.mpf('0.05')):
        fo = moments_formula(N, tt, k)[k]
        gr, _ = weyl_grid(N, [float(tt)], k, 2 * (k + N - 1) + 3)
        res.append((float(fo[0]), float(fo[1]), gr[(0, k)]))
    (tr0, fu0, g0), (tr1, fu1, g1) = res
    check('N=%d,k=%d,t*=%.6f: Haar value == FULL Hall norm (equality, not strict)' % (N, k, float(tstar)),
          abs(fu0 - g0) < 1e-10 * g0 and abs(tr0 - g0) < 1e-10 * g0,
          'grid %.12f, full %.12f; at t*+0.05: grid %.6f, full %.6f' % (g0, fu0, g1, fu1))
# Kaluza: for t>=0 the (1^k) component (-1)^{k-1} k [y^k] log sum Q^{j^2} y^j is > 0
okk = True
for Qv in (mp.mpf('1.0001'), mp.mpf('1.1'), mp.mpf(3)):
    A = [Qv ** (j * j) for j in range(41)]
    L = [mp.mpf(0)] * 41
    for kk in range(1, 41):
        L[kk] = A[kk] - mp.fsum(j * L[j] * A[kk - j] for j in range(1, kk)) / kk
    okk &= all(L[kk] > 0 for kk in range(1, 41))
check('t>=0: [y^k] log sum Q^{j^2} y^j > 0 for k<=40 (so strict inequality holds for t>=0)', okk)

# ============================================================ D. scaling limit
print('=== D. Theorem 4: K_tau(u) ===')


def gl01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (x + 1) / 2, w / 2


def simplex(r, u, n):
    if r == 1:
        return np.array([[u]]), np.array([1.0])
    x, w = gl01(n)
    X = [gg.ravel() for gg in np.meshgrid(*([x] * (r - 1)), indexing='ij')]
    W = np.prod([gg.ravel() for gg in np.meshgrid(*([w] * (r - 1)), indexing='ij')], axis=0)
    A = np.empty((len(W), r))
    rem = np.full(len(W), float(u))
    jac = np.full(len(W), float(u) ** (r - 1))
    for j in range(r - 1):
        A[:, j] = rem * X[j]
        if j < r - 2:
            jac *= (1 - X[j]) ** (r - 2 - j)
        rem = rem * (1 - X[j])
    A[:, r - 1] = rem
    return A, W * jac


A_, W_ = simplex(4, 0.7, 12)
check('simplex quadrature: vol = u^{r-1}/(r-1)!, int prod a = u^{2r-1}/(2r-1)! (r=4, d sigma = projected Lebesgue)',
      abs(W_.sum() - 0.7 ** 3 / 6) < 1e-14 and abs((W_ * A_.prod(1)).sum() - 0.7 ** 7 / math.factorial(7)) < 1e-15)


def Phi_cont(A, tau):
    r = A.shape[1]
    if r == 1:
        return np.ones(len(A))
    f = {(s, s2): np.expm1(2 * tau * A[:, s] * A[:, s2]) for s, s2 in itertools.combinations(range(r), 2)}
    tot = np.zeros(len(A))
    for G in conn_graphs(r):
        p = np.ones(len(A))
        for ed in G:
            p = p * f[ed]
        tot += p
    return tot


def Phi_cont_rec(A, tau):
    r = A.shape[1]
    full = (1 << r) - 1
    Fs = {}
    for S in range(1, full + 1):
        cols = [s for s in range(r) if S >> s & 1]
        sa = A[:, cols].sum(1)
        s2 = (A[:, cols] ** 2).sum(1)
        Fs[S] = np.exp(tau * (sa * sa - s2))
    Ph = {}
    for S in range(1, full + 1):
        low = S & -S
        v = Fs[S].copy()
        T = (S - 1) & S
        while T:
            if T & low:
                v -= Ph[T] * Fs[S ^ T]
            T = (T - 1) & S
        Ph[S] = v
    return Ph[full]


def K_terms(tau, u, rmax, nlist):
    out = []
    for r in range(1, rmax + 1):
        n = nlist[r]
        A, W = simplex(r, u, n)
        I = 0.0
        for c0 in range(0, len(W), 40000):
            Ac, Wc = A[c0:c0 + 40000], W[c0:c0 + 40000]
            Ph = Phi_cont(Ac, tau) if r <= 5 else Phi_cont_rec(Ac, tau)
            h = np.prod(np.exp(2 * tau * Ac ** 2) / Ac, axis=1) * Ph ** 2
            I += float(np.dot(Wc, h))
        out.append(u * u * math.exp(-2 * tau * u) * I / math.factorial(r))
    return out


NL = {1: 1, 2: 40, 3: 40, 4: 30, 5: 16, 6: 11, 7: 8, 8: 6}
NL2 = {1: 1, 2: 30, 3: 30, 4: 22, 5: 13, 6: 9, 7: 6, 8: 5}
tq = time.time()
Kt = K_terms(2.0, 0.7, 7, NL)
Kt2 = K_terms(2.0, 0.7, 7, NL2)
print('  tau=2, u=0.7 terms T_r:', ' '.join('%.12e' % x for x in Kt), '(%.1f s)' % (time.time() - tq))
print('  coarser quadrature diffs:', ' '.join('%.1e' % abs(a - b) for a, b in zip(Kt, Kt2)))
check('quadrature converged for r<=7 (abs 1e-12)', max(abs(a - b) for a, b in zip(Kt, Kt2)) < 1e-12)
check('r=1 term = u e^{-2 tau u (1-u)}', abs(Kt[0] - 0.7 * math.exp(-4 * 0.7 * 0.3)) < 1e-15)
A3, _ = simplex(3, 0.7, 7)
check('Phi_r: recursion == enumeration on nodes (r=3)', np.allclose(Phi_cont(A3, 2.0), Phi_cont_rec(A3, 2.0), rtol=1e-12))
A5, _ = simplex(5, 0.7, 5)
check('Phi_r: recursion == enumeration on nodes (r=5)', np.allclose(Phi_cont(A5, 2.0), Phi_cont_rec(A5, 2.0), rtol=1e-10))
K4 = sum(Kt[:4])
K6 = sum(Kt)
print('  K_2(0.7): r<=4: %.12f   r<=7: %.12f   (T_5 = %.2e, T_6 = %.2e, T_7 = %.2e)' % (K4, K6, Kt[4], Kt[5], Kt[6]))
for r in range(2, 8):
    tau, u = 2.0, 0.7
    B = (2 * tau) ** (2 * r - 2) * u ** (4 * r - 3) / (math.factorial(r) * math.factorial(2 * r - 1))
    lo, hi = math.exp(-2 * tau * u) * B, math.exp(6 * tau * u * u) * B
    check('two-sided term bound (4.6) r=%d' % r, lo <= Kt[r - 1] <= hi, '%.2e <= %.2e <= %.2e' % (lo, Kt[r - 1], hi))


# ---- finite N, exact formula, graph-free: Hall norm of q_k from the log recursion in the p-basis
def mulL(A, B):
    C = {}
    for a, x in A.items():
        for b, y in B.items():
            key = tuple(sorted(a + b, reverse=True))
            C[key] = C.get(key, 0) + x * y
    return C


def finiteN(N, k, tau, use_mp):
    t = (mp.mpf(tau) if use_mp else tau) / N ** 2
    Qv = mp.exp(t) if use_mp else math.exp(t)
    one = mp.mpf(1) if use_mp else 1.0
    g = [None] + [{lam: one * (-1) ** (j - len(lam)) / zee(lam) * Qv ** (j * j) for lam in partitions(j)}
                  for j in range(1, k + 1)]
    l = [None] * (k + 1)
    for kk in range(1, k + 1):
        acc = dict(g[kk])
        for j in range(1, kk):
            prod_ = mulL(l[j], g[kk - j])
            for key, v in prod_.items():
                acc[key] = acc.get(key, 0) - j * v / kk
        l[kk] = acc
    pref = (mp.exp(-2 * t * N * k) if use_mp else math.exp(-2 * t * N * k)) / N
    byr = Counter()
    for lam, v in l[k].items():
        c = (-1) ** (k - 1) * k * v
        byr[len(lam)] += pref * c * c * zee(lam)
    tot = mp.fsum(byr.values()) if use_mp else math.fsum(byr.values())
    return tot, byr


# consistency of the finite-N route with Theorem 1 coefficients at k=7
mp.mp.dps = 30
v7, _ = finiteN(10, 7, 2.0, True)
tN = mp.mpf(2) / 100
v7b = mp.exp(-2 * tN * 70) * mp.fsum(evalQ(CLAM[l], mp.exp(tN)) ** 2 * zee(l) for l in partitions(7)) / 10
check('finite-N value: log-recursion route == graph/CLAM route (N=10,k=7)', abs(v7 - v7b) < mp.mpf('1e-25'))


def limit_table(tau, u, Ns, Kref, label):
    vals = {}
    for N in Ns:
        k = round(u * N)
        assert abs(k - u * N) < 1e-9
        tq = time.time()
        mp.mp.dps = 40
        v, byr = finiteN(N, k, tau, True)
        if N == 40:
            vf, byrf = finiteN(N, k, tau, False)
            rmax = max(byr)
            check('%s: double vs 40-digit total at N=40' % label, abs(float(v) - vf) < 1e-12 * abs(vf),
                  '%.1e; but r=%d term: double %.1e vs 40-digit %.1e (rounding noise)' % (abs(float(v) - vf), rmax, float(byrf[rmax]), float(byr[rmax])))
        vals[N] = float(v)
        # domination: r-th term <= M_r (tau>=0) or the |tau| version with e^{2|tau|}
        at = abs(tau)
        okb = True
        for r, Tr in byr.items():
            Mr = math.exp(6 * at) * 2 ** (2 * r - 1) * (2 * at) ** (2 * r - 2) / (math.factorial(r) * math.factorial(2 * r - 1))
            if tau < 0:
                Mr *= math.exp(2 * at)
            okb &= float(Tr) <= Mr * (1 + 1e-12)
        check('%s N=%d: every r-term below the uniform bound M_r' % (label, N), okb)
        print('  %s N=%d k=%d: E|p_k|^2/N = %.12f  rel.err vs K %.2e  (%.1f s)' % (label, N, k, vals[N], vals[N] / Kref - 1, time.time() - tq))
    Ns_ = sorted(vals)
    for a, b in zip(Ns_, Ns_[1:]):
        R = (b * b * vals[b] - a * a * vals[a]) / (b * b - a * a)
        print('    Richardson 1/N^2 (%d,%d): %.12f  rel.err %.2e' % (a, b, R, R / Kref - 1))
    a, b = Ns_[-2], Ns_[-1]
    R = (b * b * vals[b] - a * a * vals[a]) / (b * b - a * a)
    err_raw = abs(vals[b] / Kref - 1)
    return R, err_raw, vals


R, er, vals = limit_table(2.0, 0.7, [10, 20, 30, 40], K6, 'tau=2,u=0.7')
check('tau=2,u=0.7: finite N -> K (raw error at N=40 < 2e-3, Richardson < 2e-5)', er < 2e-3 and abs(R / K6 - 1) < 2e-5,
      'raw %.2e, Richardson %.2e' % (er, abs(R / K6 - 1)))
e1 = abs(vals[20] / K6 - 1)
e2 = abs(vals[40] / K6 - 1)
print('  error ratio N=20 -> 40: %.3f (1/N^2 predicts 4)' % (e1 / e2))

# negative tau
Ktn = K_terms(-1.5, 0.6, 7, NL)
Ktn2 = K_terms(-1.5, 0.6, 7, NL2)
check('tau=-1.5 quadrature converged', max(abs(a - b) for a, b in zip(Ktn, Ktn2)) < 1e-12)
Kn = sum(Ktn)
print('  K_{-1.5}(0.6) = %.12f; terms' % Kn, ' '.join('%.3e' % x for x in Ktn))
R, er, _ = limit_table(-1.5, 0.6, [10, 20, 30, 40], Kn, 'tau=-1.5,u=0.6')
check('tau=-1.5,u=0.6: finite N -> K (Richardson < 2e-5)', abs(R / Kn - 1) < 2e-5, 'raw %.2e, Richardson %.2e' % (er, abs(R / Kn - 1)))

# author's table value tau=2, u=1
KaT = K_terms(2.0, 1.0, 8, NL)
KaT2 = K_terms(2.0, 1.0, 8, NL2)
Ka = sum(KaT)
print('  tau=2,u=1 terms r=6,7,8: %.3e %.3e %.3e; coarse diffs %s' % (KaT[5], KaT[6], KaT[7], ' '.join('%.0e' % abs(a - b) for a, b in zip(KaT, KaT2))))
check('author table K_2(1) = 1.76699977 (r<=8)', abs(Ka - 1.76699977) < 1e-8, '%.10f (r<=6 alone: %.10f)' % (Ka, sum(KaT[:6])))

# small-tau expansion, second order (claimed) and third order (derived here)
for u in (0.7, 0.35, 1.0):
    C2 = 2 * u ** 3 * (1 - u) ** 2 + u ** 5 / 3
    C3 = -4 / 3 * u ** 4 * (1 - u) ** 3 + 8 / 15 * u ** 7 - 2 / 3 * u ** 6
    est2, est3 = [], []
    for tau in (0.02, 0.01):
        Kp = sum(K_terms(tau, u, 4, NL))
        Km = sum(K_terms(-tau, u, 4, NL))
        est2.append(((Kp + Km) / 2 - u) / tau ** 2)
        est3.append(((Kp - Km) / 2 + 2 * tau * u * u * (1 - u)) / tau ** 3)
    c2 = (4 * est2[1] - est2[0]) / 3
    c3 = (4 * est3[1] - est3[0]) / 3
    check('u=%g: tau^2 coefficient = 2u^3(1-u)^2 + u^5/3' % u, abs(c2 - C2) < 1e-7, '%.9f vs %.9f' % (c2, C2))
    check('u=%g: tau^3 coefficient = -4/3 u^4(1-u)^3 + 8/15 u^7 - 2/3 u^6' % u, abs(c3 - C3) < 1e-5, '%.7f vs %.7f' % (c3, C3))

print('time %.1f s' % (time.time() - T0))
print('%d checks, %d failed' % (NCHK[0], NCHK[1]))
