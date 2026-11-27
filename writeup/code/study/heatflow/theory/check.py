#!/usr/bin/env python3
"""check.py -- recompute and assert every number in proof.md (flowed CUE power sums).

Run:  (ulimit -v 1500000; timeout 600 nice /home/bas/praime/.venv/bin/python3 -I check.py)

Notation (as in proof.md):  E_t(x) = sum_j exp(-t j (N-j)) e_j x^j = prod_i (1 + zeta_i x),
p_k(t) = sum_i zeta_i^k,  Q = e^t,  F_t(y) = sum_j Q^{j^2} e_j y^j,  q_k = e^{tNk} p_k(t).
Claim 1:  q_k = sum_{lam |- k} c_lam p_lam,
          c_lam = (-1)^{r-1} (k / z_lam) Q^{sum i_s^2} Phi_lam,
          Phi_lam = sum_{connected G on [r]} prod_{edges} (Q^{2 i_s i_s'} - 1).
Claim 2:  E|p_k(t)|^2 = e^{-2tNk} H_k,  H_k = sum_lam c_lam^2 z_lam   (k <= N);
          for all k:  E|p_k(t)|^2 = e^{-2tNk} sum_{nu |- k, l(nu) <= N} (sum_lam c_lam chi^nu_lam)^2.
Claim 3:  E|p_k|^2 / N -> K_tau(u)  (t = tau/N^2, k = uN).
"""
import itertools, math, sys, time
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import numpy as np
import sympy as sp

T0 = time.time()
NCHK = 0
NFAIL = 0


def check(name, ok, detail=""):
    global NCHK, NFAIL
    NCHK += 1
    if not ok:
        NFAIL += 1
    print(f"[{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    sys.stdout.flush()


# ---------------------------------------------------------------- partitions
def partitions(n, m=None):
    if m is None:
        m = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, m), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def zee(lam):
    z = 1
    for i, m in Counter(lam).items():
        z *= i ** m * math.factorial(m)
    return z


def merge(a, b):
    return tuple(sorted(a + b, reverse=True))


# ---------------------------------------------------------------- Route A: log recursion in the p-basis
def routeA(kmax, coef, scal):
    """Return dict k -> {lam: c_lam}, from q_k = (-1)^{k-1} k [y^k] log(sum_j coef(j) e_j y^j).
    coef(j) is the multiplier (Q^{j^2} or exp(t j^2)); scal(x) converts a Fraction to the ring."""
    g = {}
    for j in range(1, kmax + 1):
        cj = coef(j)
        g[j] = {lam: cj * scal(Fraction((-1) ** (j - len(lam)), zee(lam))) for lam in partitions(j)}
    L = {}
    for n in range(1, kmax + 1):
        acc = {lam: v * scal(Fraction(n)) for lam, v in g[n].items()}
        for m in range(1, n):
            for la, va in L[m].items():
                for lb, vb in g[n - m].items():
                    key = merge(la, lb)
                    term = va * vb * scal(Fraction(-m))
                    acc[key] = acc[key] + term if key in acc else term
        L[n] = {lam: v * scal(Fraction(1, n)) for lam, v in acc.items()}
    return {k: {lam: v * scal(Fraction((-1) ** (k - 1) * k)) for lam, v in L[k].items()}
            for k in range(1, kmax + 1)}


# ---------------------------------------------------------------- Route B: literal connected graphs
_CONN = {}


def connected_graphs(r):
    """List of edge lists (pairs) of all connected graphs on vertex set {0..r-1}."""
    if r in _CONN:
        return _CONN[r]
    edges = list(itertools.combinations(range(r), 2))
    out = []
    for mask in range(1 << len(edges)):
        par = list(range(r))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        es = [edges[i] for i in range(len(edges)) if mask >> i & 1]
        for a, b in es:
            par[find(a)] = find(b)
        if len({find(x) for x in range(r)}) == 1:
            out.append(es)
    _CONN[r] = out
    return out


Qs = sp.Symbol('Q')


def P(e):
    return sp.Poly(e, Qs, domain='QQ')


def Phi_literal(lam):
    r = len(lam)
    cnt = Counter()
    for es in connected_graphs(r):
        cnt[tuple(sorted(2 * lam[a] * lam[b] for a, b in es))] += 1
    tot = P(0)
    for key, c in cnt.items():
        term = P(c)
        for w in key:
            term = term * P(Qs ** w - 1)
        tot = tot + term
    return tot


def Phi_recursive(parts, f, A):
    """Connected sum over graphs on len(parts) labelled vertices with edge weight f, via
    A(S) = sum_{set partitions} prod conn(blocks)  (A = all-graph sum = prod (1+f))."""
    r = len(parts)
    conn = {}
    for S in range(1, 1 << r):
        low = S & -S
        rest = S ^ low
        val = A(S)
        sub = rest
        while True:          # enumerate proper subsets T = low|sub of S
            if sub != rest:
                T = low | sub
                val = val - conn[T] * A(S ^ T)
            if sub == 0:
                break
            sub = (sub - 1) & rest
        conn[S] = val
    return conn[(1 << r) - 1]


def c_B(lam, Phi):
    r, k = len(lam), sum(lam)
    return P(sp.Rational((-1) ** (r - 1) * k, zee(lam)) * Qs ** sum(i * i for i in lam)) * Phi


# ---------------------------------------------------------------- characters (Murnaghan-Nakayama)
@lru_cache(maxsize=None)
def chi(nu, lam):
    if not lam:
        return 1 if sum(nu) == 0 else 0
    h, rest = lam[0], lam[1:]
    Lh = len(nu)
    beta = [nu[i] + (Lh - 1 - i) for i in range(Lh)]
    bs = set(beta)
    tot = 0
    for b in beta:
        if b - h >= 0 and (b - h) not in bs:
            sign = (-1) ** sum(1 for c in beta if b - h < c < b)
            nb = sorted([c for c in beta if c != b] + [b - h], reverse=True)
            nn = tuple(x for x in (nb[i] - (Lh - 1 - i) for i in range(Lh)) if x > 0)
            tot += sign * chi(nn, rest)
    return tot


def truncated_norm(c, k, N):
    """sum over nu |- k with l(nu) <= N of (sum_lam c_lam chi^nu(lam))^2."""
    tot = 0
    for nu in partitions(k):
        if len(nu) <= N:
            s = sum(v * chi(nu, lam) for lam, v in c.items())
            tot += s * s
    return tot


def hall_norm(c):
    return sum(v * v * zee(lam) for lam, v in c.items())


print("=== 0. implementation sanity: MN characters, column orthogonality, k <= 8")
ok = True
for k in range(1, 9):
    pk = list(partitions(k))
    for la in pk:
        for mu in pk:
            s = sum(chi(nu, la) * chi(nu, mu) for nu in pk)
            ok &= (s == (zee(la) if la == mu else 0))
check("sum_nu chi^nu(lam) chi^nu(mu) = delta z_lam for all lam, mu |- k <= 8", ok)

# ================================================================ Claim 1, exact
print("\n=== 1. Claim 1: c_lam from the log recursion (route A) = cluster formula (route B), exact in Q = e^t")
KMAX_EX = 8
cA = routeA(KMAX_EX, lambda j: P(Qs ** (j * j)), lambda fr: P(sp.Rational(fr.numerator, fr.denominator)))
nlit = nrec = 0
ok_lit = ok_rec = True
for k in range(1, KMAX_EX + 1):
    for lam in partitions(k):
        r = len(lam)
        if r <= 6:
            Phi = Phi_literal(lam)
            ok_lit &= (c_B(lam, Phi) == cA[k][lam])
            nlit += 1
        else:
            Phi = Phi_recursive(lam, None, lambda S, lam=lam: P(Qs ** (2 * sum(
                lam[a] * lam[b] for a in range(len(lam)) for b in range(a + 1, len(lam))
                if S >> a & 1 and S >> b & 1))))
            ok_rec &= (c_B(lam, Phi) == cA[k][lam])
            nrec += 1
check(f"route A = route B (literal enumeration of connected graphs, r <= 6), all {nlit} partitions lam |- k <= {KMAX_EX} with l(lam) <= 6", ok_lit)
check(f"route A = route B (connected sum by set-partition recursion, r = 7, 8), {nrec} partitions", ok_rec)
# literal vs recursive Phi agree (tests the recursion itself), all lam |- k <= 7 with r <= 6
ok = True
for k in range(1, 8):
    for lam in partitions(k):
        if len(lam) <= 6:
            Phr = Phi_recursive(lam, None, lambda S, lam=lam: P(Qs ** (2 * sum(
                lam[a] * lam[b] for a in range(len(lam)) for b in range(a + 1, len(lam))
                if S >> a & 1 and S >> b & 1))))
            ok &= (Phr == Phi_literal(lam))
check("literal connected-graph sum = set-partition recursion, all lam |- k <= 7, l(lam) <= 6", ok)
# number of connected labelled graphs (OEIS A001187: 1, 1, 4, 38, 728, 26704)
check("connected labelled graphs on r = 1..6 vertices: 1, 1, 4, 38, 728, 26704",
      [len(connected_graphs(r)) for r in range(1, 7)] == [1, 1, 4, 38, 728, 26704])
# t = 0 (Q = 1): q_k = p_k
ok = all(cA[k][lam].eval(1) == (1 if lam == (k,) else 0) for k in range(1, KMAX_EX + 1) for lam in partitions(k))
check("at t = 0: c_lam = delta_{lam,(k)}", ok)
# explicit low cases
check("k = 2: q_2 = Q^4 p_2 - (Q^4 - Q^2) p_1^2",
      cA[2][(2,)] == P(Qs ** 4) and cA[2][(1, 1)] == P(-(Qs ** 4 - Qs ** 2)))
ok = all(cA[k][(k,)] == P(Qs ** (k * k)) for k in range(1, KMAX_EX + 1))
check("c_(k) = Q^{k^2} = e^{t k^2} exactly, k <= 8", ok)
# first-order term (second route, via D^2 E / E):  q_k = (1 + t k^2) p_k - t k sum_{a=1}^{k-1} p_a p_{k-a} + O(t^2)
ok = True
for k in range(2, KMAX_EX + 1):
    for lam in partitions(k):
        cp = cA[k][lam]
        d1 = sum(coef * d for (d,), coef in zip(cp.monoms(), cp.coeffs()))   # d/dt at t=0 of sum h_d e^{dt}
        if lam == (k,):
            want = k * k
        elif len(lam) == 2:
            a, b = lam
            want = -k * (2 if a != b else 1)
        else:
            want = 0
        ok &= (d1 == want)
check("d/dt c_lam at t = 0 equals the D^2E/E route: k^2 for (k), -k*(#ordered (a,b)) for 2 parts, 0 for r >= 3", ok)

# float implementation of route A (used below) = exact route A at Q = e^t, k <= 8
cf = routeA(8, lambda j: math.exp(0.05 * j * j), lambda fr: float(fr))
errf = max(abs(cf[k][lam] - float(cA[k][lam].eval(sp.exp(sp.Rational(1, 20))).evalf(30))) / max(1.0, abs(cf[k][lam]))
           for k in range(1, 9) for lam in partitions(k))
check("float route A = exact route A at t = 0.05, all lam |- k <= 8", errf < 1e-12, f"max err {errf:.1e}")

# ================================================================ Claim 1 numerically at finite N, including k > N
print("\n=== 1b. Claim 1 at N = 6 for a random U, k = 1..10 (k > N included), several t")
from scipy.stats import unitary_group
rng = np.random.default_rng(20261127)
cAf_cache = {}


def cA_float(kmax, t):
    key = (kmax, t)
    if key not in cAf_cache:
        cAf_cache[key] = routeA(kmax, lambda j: math.exp(t * j * j), lambda fr: float(fr))
    return cAf_cache[key]


def newton_from_coeffs(a, kmax):
    """a[0]=1, a[j] = e_j(zeta) (array over samples); return p_1..p_kmax of zeta."""
    Nn = len(a) - 1
    p = [None] * (kmax + 1)
    for k in range(1, kmax + 1):
        s = (-1) ** (k - 1) * k * (a[k] if k <= Nn else 0)
        for j in range(1, k):
            if j <= Nn:
                s = s + (-1) ** (j - 1) * a[j] * p[k - j]
        p[k] = s
    return p


def elem(z):
    """z shape (M, N) -> list e_0..e_N arrays of shape (M,)."""
    M, Nn = z.shape
    e = [np.ones(M, complex)] + [np.zeros(M, complex) for _ in range(Nn)]
    for i in range(Nn):
        for j in range(i + 1, 0, -1):
            e[j] = e[j] + z[:, i] * e[j - 1]
    return e


U = unitary_group.rvs(6, random_state=rng)
z1 = np.linalg.eigvals(U)[None, :]
N1 = 6
e1 = elem(z1)
worst = 0.0
for t in (0.03, -0.02, 0.2):
    a = [e1[j] * math.exp(-t * j * (N1 - j)) for j in range(N1 + 1)]
    pt = newton_from_coeffs(a, 10)
    c = cA_float(10, t)
    for k in range(1, 11):
        val = sum(v * np.prod([np.sum(z1[0] ** i) for i in lam]) for lam, v in c[k].items())
        val *= math.exp(-t * N1 * k)
        worst = max(worst, abs(val - pt[k][0]) / max(1, abs(pt[k][0])))
check("p_k(t) (Newton from E_t) = e^{-tNk} sum c_lam p_lam, N = 6, k = 1..10, t in {0.03,-0.02,0.2}",
      worst < 1e-9, f"max rel err {worst:.1e}")
# Newton power sums = power sums of the actual roots
t = 0.03
a = [e1[j][0] * math.exp(-t * j * (N1 - j)) for j in range(N1 + 1)]
xr = np.roots(a[::-1])
zeta = -1 / xr
pt = newton_from_coeffs([np.array([v]) for v in a], 8)
err = max(abs(np.sum(zeta ** k) - pt[k][0]) for k in range(1, 9))
check("Newton power sums = power sums of -1/(roots of E_t), N = 6, t = 0.03, k <= 8", err < 1e-9, f"{err:.1e}")

# ================================================================ Claim 2, exact via Weyl integration
print("\n=== 2. Claim 2 exact: Weyl integration (constant term) at rational Q, q_k from Newton (no c_lam used)")


def poly_mul(A, B):
    C = {}
    for ea, va in A.items():
        for eb, vb in B.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            C[e] = C.get(e, 0) + va * vb
    return {e: v for e, v in C.items() if v != 0}


def poly_add(A, B, s=1):
    C = dict(A)
    for e, v in B.items():
        C[e] = C.get(e, 0) + s * v
    return {e: v for e, v in C.items() if v != 0}


def weyl_expect_abs2(q, Nn):
    V = {tuple([0] * Nn): Fraction(1)}
    for i in range(Nn):
        for j in range(Nn):
            if i != j:
                e = [0] * Nn
                e[i] += 1
                e[j] -= 1
                V = poly_mul(V, {tuple([0] * Nn): Fraction(1), tuple(e): Fraction(-1)})
    qV = poly_mul(q, V)
    return sum(v * qV.get(e, 0) for e, v in q.items()) / math.factorial(Nn)


def q_newton_exact(Nn, Qv, kmax):
    ej = []
    for j in range(Nn + 1):
        d = {}
        for S in itertools.combinations(range(Nn), j):
            e = tuple(1 if i in S else 0 for i in range(Nn))
            d[e] = Fraction(1)
        ej.append(d)
    a = [{e: v * Qv ** (j * j) for e, v in ej[j].items()} for j in range(Nn + 1)]
    q = [None] * (kmax + 1)
    for k in range(1, kmax + 1):
        s = {e: (-1) ** (k - 1) * k * v for e, v in a[k].items()} if k <= Nn else {}
        for j in range(1, min(k - 1, Nn) + 1):
            s = poly_add(s, poly_mul(a[j], q[k - j]), (-1) ** (j - 1))
        q[k] = s
    return q


weyl_rows = []
for Qv in (Fraction(11, 10), Fraction(4, 5)):
    QR = sp.Rational(Qv.numerator, Qv.denominator)
    cAQ = {}
    for k in range(1, 7):
        cAQ[k] = {}
        for lam, v in cA[k].items():
            val = sp.Rational(v.eval(QR))
            cAQ[k][lam] = Fraction(int(val.p), int(val.q))
    for Nn, kmax in ((2, 4), (3, 5), (4, 6)):
        qs = q_newton_exact(Nn, Qv, kmax)
        for k in range(1, kmax + 1):
            W = weyl_expect_abs2(qs[k], Nn)
            Hk = hall_norm(cAQ[k])
            Tk = truncated_norm(cAQ[k], k, Nn)
            weyl_rows.append((Qv, Nn, k, W, Hk, Tk))
ok_le = all(W == Hk for Qv, Nn, k, W, Hk, Tk in weyl_rows if k <= Nn)
ok_tr = all(W == Tk for Qv, Nn, k, W, Hk, Tk in weyl_rows)
ok_gt = all(W != Hk for Qv, Nn, k, W, Hk, Tk in weyl_rows if k > Nn)
n_le = sum(1 for r_ in weyl_rows if r_[2] <= r_[1])
check(f"k <= N: E|q_k|^2 (Weyl, exact) = sum c_lam^2 z_lam exactly, {n_le} cases (Q = 11/10, 4/5; N = 2,3,4)", ok_le)
check(f"all k: E|q_k|^2 (Weyl) = truncated Schur norm over l(nu) <= N exactly, {len(weyl_rows)} cases", ok_tr)
check("k > N: full Hall norm is NOT the expectation (strict inequality in every case)", ok_gt)
for Qv, Nn, k, W, Hk, Tk in weyl_rows:
    if Qv == Fraction(11, 10) and k >= Nn:
        print(f"    Q=11/10 N={Nn} k={k}: Weyl {float(W):.6f}  Hall {float(Hk):.6f}  truncated {float(Tk):.6f}")

# ================================================================ circle property and zero dynamics
print("\n=== 3. roots of E_t on |x| = 1 for t >= 0 (Lee-Yang); off the circle for t < 0; zero dynamics")
Ns = 8
Us = unitary_group.rvs(Ns, size=2000, random_state=rng)
zs = np.linalg.eigvals(Us)
es = elem(zs)
for t in (0.001, 0.03, 0.3):
    dev = 0.0
    for m in range(zs.shape[0]):
        a = [es[j][m] * math.exp(-t * j * (Ns - j)) for j in range(Ns + 1)]
        dev = max(dev, np.max(np.abs(np.abs(np.roots(a[::-1])) - 1)))
    check(f"t = {t}: all roots of E_t on the unit circle, 2000 Haar U (N = 8)", dev < 1e-7, f"max ||x|-1| = {dev:.1e}")
t = -0.03
noff = 0
for m in range(zs.shape[0]):
    a = [es[j][m] * math.exp(-t * j * (Ns - j)) for j in range(Ns + 1)]
    if np.max(np.abs(np.abs(np.roots(a[::-1])) - 1)) > 1e-6:
        noff += 1
aI = [math.comb(Ns, j) * math.exp(-t * j * (Ns - j)) for j in range(Ns + 1)]
devI = np.max(np.abs(np.abs(np.roots(aI[::-1])) - 1))
check("t = -0.03, U = I (N = 8): some root off the circle", devI > 0.05, f"max ||x|-1| = {devI:.3f}")
print(f"    t = -0.03: {noff} of 2000 Haar U have a root off the circle (|.|-1 > 1e-6)")
check("t = -0.03: a positive fraction of Haar U have roots off the circle", noff > 0, f"{noff}/2000")
# zero dynamics d theta_i/dt = sum_{j != i} cot((theta_i - theta_j)/2)
m = 7
h = 1e-6
t = 0.02


def roots_at(t):
    a = [es[j][m] * math.exp(-t * j * (Ns - j)) for j in range(Ns + 1)]
    return np.roots(a[::-1])


x0, xp, xm = roots_at(t), roots_at(t + h), roots_at(t - h)
th0 = np.angle(x0)
num = np.array([np.angle(xp[np.argmin(abs(xp - x))] / xm[np.argmin(abs(xm - x))]) / (2 * h) for x in x0])
pred = np.array([sum(1 / math.tan((th0[i] - th0[j]) / 2) for j in range(Ns) if j != i) for i in range(Ns)])
err = np.max(np.abs(num - pred)) / np.max(np.abs(pred))
check("zero dynamics: d theta_i/dt = sum_{j!=i} cot((theta_i-theta_j)/2) (finite differences, N = 8, t = 0.02)",
      err < 1e-5, f"rel err {err:.1e}")

# consequence of circle + Claim 2 at N = k:  H_k(t) <= k^2 e^{2 t k^2} for t >= 0 (not an identity)
ok = True
ratios = []
for t in (0.01, 0.1, 0.3):
    c = cA_float(12, t)
    for k in range(1, 13):
        rat = hall_norm(c[k]) / (k * k * math.exp(2 * t * k * k))
        ratios.append((t, k, rat))
        ok &= rat <= 1 + 1e-12
check("H_k(t) <= k^2 e^{2tk^2}, t in {0.01,0.1,0.3}, k <= 12", ok,
      "ratio at t=0.3,k=6: %.6f; t=0.3,k=12: %.6f" % ([r_[2] for r_ in ratios if r_[0] == 0.3 and r_[1] == 6][0],
                                                        [r_[2] for r_ in ratios if r_[0] == 0.3 and r_[1] == 12][0]))
# exact version at rational Q >= 1 for k <= 6: H_k(Q) <= k^2 Q^{2k^2}
ok = True
for Qv in (sp.Rational(11, 10), sp.Rational(3, 2), sp.Rational(5, 1)):
    for k in range(1, 7):
        Hk = sum(cA[k][lam].eval(Qv) ** 2 * zee(lam) for lam in partitions(k))
        ok &= Hk <= k * k * Qv ** (2 * k * k)
check("exact: H_k(Q) <= k^2 Q^{2k^2} at Q = 11/10, 3/2, 5 (k <= 6)", ok)

# ================================================================ small-t expansion of H_k, exact
print("\n=== 4. H_k(t) = k e^{2tk^2} + t^2 k^2 (k^3-k)/3 + O(t^3), exact for k <= 8")
ok = True
for k in range(1, KMAX_EX + 1):
    H = P(0)
    for lam in partitions(k):
        H = H + cA[k][lam] ** 2 * P(zee(lam))
    mon = [d for (d,) in H.monoms()]
    co = H.coeffs()
    tay = [sum(cc * sp.Integer(d) ** mm for d, cc in zip(mon, co)) / sp.factorial(mm) for mm in range(3)]
    ok &= tay[0] == k and tay[1] == 2 * k ** 3 and tay[2] == 2 * k ** 5 + sp.Rational(k * k * (k ** 3 - k), 3)
check("Taylor coefficients of H_k(e^t): [k, 2k^3, 2k^5 + k^2(k^3-k)/3], k = 1..8", ok)

# ================================================================ Claim 3: K_tau(u)
print("\n=== 5. Claim 3: the limit K_tau(u)")


def phi_vec(a, tau):
    """a: array (r, M). Connected-graph sum with f = exp(2 tau a_s a_s') - 1, via set partitions."""
    r = a.shape[0]

    def A(S):
        idx = [s for s in range(r) if S >> s & 1]
        sub = a[idx]
        return np.exp(tau * (sub.sum(0) ** 2 - (sub ** 2).sum(0)))
    Acache = {}

    def Ac(S):
        if S not in Acache:
            Acache[S] = A(S)
        return Acache[S]
    conn = {}
    for S in range(1, 1 << r):
        low = S & -S
        rest = S ^ low
        val = Ac(S).copy()
        sub = rest
        while True:
            if sub != rest:
                T = low | sub
                val -= conn[T] * Ac(S ^ T)
            if sub == 0:
                break
            sub = (sub - 1) & rest
        conn[S] = val
    return conn[(1 << r) - 1]


def K_terms(tau, u, rmax=7, n=None):
    """T_r = (1/r!) int_{simplex sum a = u} prod(e^{2 tau a^2}/a) Phi_r^2 da_1..da_{r-1}."""
    terms = [math.exp(2 * tau * u * u) / u]
    for r in range(2, rmax + 1):
        nn = n if n else {2: 40, 3: 30, 4: 18, 5: 11, 6: 8, 7: 6}[r]
        x, w = np.polynomial.legendre.leggauss(nn)
        x = (x + 1) / 2
        w = w / 2
        grids = np.meshgrid(*([x] * (r - 1)), indexing='ij')
        wg = np.meshgrid(*([w] * (r - 1)), indexing='ij')
        X = [g.ravel() for g in grids]
        W = np.prod([g.ravel() for g in wg], axis=0)
        rem = np.full(X[0].shape, u)
        a = []
        J = np.ones_like(rem)
        for s in range(r - 1):
            a.append(rem * X[s])
            J *= rem
            rem = rem * (1 - X[s])
        a.append(rem)
        a = np.array(a)
        ph = phi_vec(a, tau)
        h = np.prod(np.exp(2 * tau * a ** 2) / a, axis=0) * ph ** 2
        terms.append(float(np.sum(W * J * h)) / math.factorial(r))
    return terms


def K(tau, u, **kw):
    T = K_terms(tau, u, **kw)
    return u * u * math.exp(-2 * tau * u) * sum(T), T


# quadrature accuracy: two node counts
ok = True
for tau, u in ((1.0, 0.75), (0.5, 1.0)):
    T1 = K_terms(tau, u, rmax=4)
    T2 = K_terms(tau, u, rmax=4, n=24)
    ok &= all(abs(p - q) <= 1e-10 * max(1e-300, abs(q)) + 1e-15 for p, q in zip(T1, T2))
check("simplex quadrature stable under change of node count (r <= 4)", ok)

# tree bounds on Phi_r:  (2tau)^{r-1} u^{r-2} prod a <= Phi_r <= e^{2 tau u^2} (2tau)^{r-1} u^{r-2} prod a
ok = True
for r in range(2, 7):
    for tau in (0.3, 1.0, 3.0):
        a = rng.dirichlet(np.ones(r), size=2000).T * 0.8
        ph = phi_vec(a, tau)
        base = (2 * tau) ** (r - 1) * 0.8 ** (r - 2) * np.prod(a, axis=0)
        ok &= bool(np.all(ph >= base * (1 - 1e-9)) and np.all(ph <= math.exp(2 * tau * 0.64) * base * (1 + 1e-9)))
check("tree bounds on Phi_r at 2000 random points, r = 2..6, tau in {0.3,1,3}, u = 0.8", ok)
# weighted Cayley: sum_T prod x^{deg} = prod x (sum x)^{r-2}, and compositions: sum prod i_s = C(k+r-1, 2r-1)
ok = True
for r in range(2, 6):
    xs = rng.random(r) + 0.1
    edges = list(itertools.combinations(range(r), 2))
    tot = 0.0
    for es_ in itertools.combinations(edges, r - 1):
        par = list(range(r))

        def fd(v):
            while par[v] != v:
                v = par[v]
            return v
        good = True
        for p_, q_ in es_:
            a_, b_ = fd(p_), fd(q_)
            if a_ == b_:
                good = False
                break
            par[a_] = b_
        if good:
            deg = Counter([v for e_ in es_ for v in e_])
            tot += np.prod([xs[v] ** deg[v] for v in range(r)])
    ok &= abs(tot - np.prod(xs) * xs.sum() ** (r - 2)) < 1e-10 * tot
check("weighted Cayley formula sum_T prod x_s^{deg_T(s)} = prod x (sum x)^{r-2}, r = 2..5", ok)
ok = True
for k in range(1, 13):
    for r in range(1, k + 1):
        s = 0
        for comp in itertools.combinations(range(1, k), r - 1):
            cuts = (0,) + comp + (k,)
            s += math.prod(cuts[i + 1] - cuts[i] for i in range(r))
        ok &= s == math.comb(k + r - 1, 2 * r - 1)
check("sum over compositions of k into r parts of prod i_s = C(k+r-1, 2r-1), k <= 12", ok)

# small-tau expansion of K, route (a): from the integral
for u in (0.4, 0.75, 1.0):
    def g(tau):
        return (K(tau, u, rmax=3)[0] - u + 2 * tau * u * u * (1 - u)) / tau ** 2
    est = 2 * g(1e-3) - g(2e-3)
    want = 2 * u ** 3 * (1 - u) ** 2 + u ** 5 / 3
    check(f"tau^2 coefficient of K at u = {u}: 2u^3(1-u)^2 + u^5/3 = {want:.6f}", abs(est - want) < 1e-5,
          f"numerical {est:.6f}")

# finite N exact (Claim 2 formula, floats) vs K
print("    tau    u     K_tau(u)    terms r=1..7 (u^2 e^{-2tau u} T_r)            N: E|p_k|^2/N (exact, Claim 2)")
lim_rows = []
for tau in (0.5, 1.0, 2.0):
    for u, Nlist in ((0.5, (8, 16, 32, 48)), (0.75, (8, 16, 32)), (1.0, (6, 12, 24))):
        Kv, T = K(tau, u)
        pref = u * u * math.exp(-2 * tau * u)
        vals = []
        for Nn in Nlist:
            k = round(u * Nn)
            t = tau / Nn ** 2
            c = cA_float(k, t)[k]
            vals.append(math.exp(-2 * t * Nn * k) * hall_norm(c) / Nn)
        lim_rows.append((tau, u, Kv, [pref * x for x in T], list(zip(Nlist, vals))))
        print(f"    {tau:.1f}  {u:.2f}  {Kv:.8f}  " + " ".join(f"{pref * x:.1e}" for x in T) + "   " +
              "  ".join(f"N={Nn}: {v:.8f}" for Nn, v in zip(Nlist, vals)))
# convergence: error ~ C/N^2 (Euler-Maclaurin: integrand vanishes on the faces)
ok_conv = True
ok_rate = True
rates = []
for tau, u, Kv, T, NV in lim_rows:
    errs = [abs(v - Kv) for Nn, v in NV]
    ok_conv &= errs[-1] < 3e-3 * Kv and all(errs[i + 1] < errs[i] for i in range(len(errs) - 1))
    rate = math.log(errs[-2] / errs[-1]) / math.log(NV[-1][0] / NV[-2][0])
    rates.append(rate)
    ok_rate &= 1.7 < rate < 2.3
check("E|p_k|^2/N -> K_tau(u), errors decreasing in N, last error < 0.3% of K (9 cases; first run predicted 0.2%%, failed at tau=2,u=1)".replace("%%","%"), ok_conv,
      "max rel err at largest N %.2e" % max(abs(NV[-1][1] - Kv) / Kv for tau, u, Kv, T, NV in lim_rows))
check("observed convergence order between the two largest N in (1.7, 2.3) (O(1/N^2))", ok_rate,
      "orders " + ", ".join(f"{r_:.2f}" for r_ in rates))
ok = True
for tau, u, Kv, T, NV in lim_rows:
    (N1_, v1), (N2_, v2) = NV[-2], NV[-1]
    rich = (N2_ ** 2 * v2 - N1_ ** 2 * v1) / (N2_ ** 2 - N1_ ** 2)
    ok &= abs(rich - Kv) < 2e-5 * Kv + abs(v2 - Kv) / 20
    print(f"    tau={tau} u={u}: raw rel err at N={N2_}: {abs(v2 - Kv) / Kv:.2e}; Richardson (N={N1_},{N2_}) rel err {abs(rich - Kv) / Kv:.2e}")
check("Richardson extrapolation in 1/N^2 matches K_tau(u) (error < 1/20 of the raw error)", ok)
# truncation of the r-series: last computed term small
ok = all(T[-1] < 1e-6 * Kv for tau, u, Kv, T, NV in lim_rows)
check("r = 7 term < 1e-6 K in all 9 cases (series truncation negligible)", ok,
      "max ratio %.1e" % max(T[-1] / Kv for tau, u, Kv, T, NV in lim_rows))

# two-sided bounds on the terms of K (proof.md, 4.6), and the N-uniform domination M_r (4.5)
ok = True
for tau, u, Kv, T, NV in lim_rows:
    for r in range(2, len(T) + 1):
        base = (2 * tau) ** (2 * r - 2) * u ** (4 * r - 3) / (math.factorial(r) * math.factorial(2 * r - 1))
        ok &= math.exp(-2 * tau * u) * base * (1 - 1e-6) <= T[r - 1] <= math.exp(6 * tau * u * u) * base
check("e^{-2tau u} B_r <= T_r <= e^{6tau u^2} B_r, B_r = (2tau)^{2r-2}u^{4r-3}/(r!(2r-1)!), r = 2..7, 9 cases", ok)
# 60-digit arithmetic here: in double precision the r-th terms for large r are rounding noise
# (first run: noise 1e-33..1e-37 against true values down to 1e-69, so the check "failed")
import mpmath as mp
mp.mp.dps = 60
ok = True
worst = 0.0
okH = True
for tau, Nn, k in ((1.0, 16, 12), (2.0, 16, 16), (-1.0, 16, 12), (0.5, 24, 18)):
    t = tau / Nn ** 2
    tm = mp.mpf(tau) / Nn ** 2
    c = routeA(k, lambda j: mp.exp(tm * j * j), lambda fr: mp.mpf(fr.numerator) / fr.denominator)[k]
    Hm = hall_norm(c)
    Hf = hall_norm(cA_float(k, t)[k])
    okH &= abs(Hf - float(Hm)) < 1e-12 * float(Hm)
    ub = k / Nn
    for r in range(1, k + 1):
        term = float(mp.exp(-2 * tm * Nn * k) * sum(v * v * zee(lam) for lam, v in c.items() if len(lam) == r) / Nn)
        at = abs(tau)
        Mr = math.exp(2 * at * ub) * math.exp(6 * at * ub * ub) * 2 ** (2 * r - 1) * (2 * at) ** (2 * r - 2) * ub ** (4 * r - 3) \
            / (math.factorial(r) * math.factorial(2 * r - 1))
        ok &= term <= Mr * (1 + 1e-9)
        worst = max(worst, term / Mr)
check("finite-N r-th term <= M_r (uniform domination, 60 digits), 4 cases incl. tau < 0", ok, f"max term/M_r = {worst:.2e}")
check("double-precision H_k (used in the limit table) = 60-digit H_k to 1e-12, same 4 cases", okH)

# ================================================================ Monte Carlo
print("\n=== 6. Monte Carlo, Haar U (scipy.stats.unitary_group), E|p_k(t)|^2 vs Claim 2")
M = 200000
CH = 20000
mc_rows = []
for Nn, kmax in ((6, 8), (8, 6)):
    tl = (0.01, 0.03, -0.01)
    S1 = {(t, k): 0.0 for t in tl for k in range(1, kmax + 1)}
    S2 = dict(S1)
    for b in range(M // CH):
        Ub = unitary_group.rvs(Nn, size=CH, random_state=rng)
        zb = np.linalg.eigvals(Ub)
        eb = elem(zb)
        for t in tl:
            a = [eb[j] * math.exp(-t * j * (Nn - j)) for j in range(Nn + 1)]
            pt = newton_from_coeffs(a, kmax)
            for k in range(1, kmax + 1):
                v = np.abs(pt[k]) ** 2
                S1[(t, k)] += v.sum()
                S2[(t, k)] += (v * v).sum()
    for t in tl:
        c = cA_float(kmax, t)
        for k in range(1, kmax + 1):
            mean = S1[(t, k)] / M
            se = math.sqrt((S2[(t, k)] / M - mean ** 2) / (M - 1))
            pref = math.exp(-2 * t * Nn * k)
            full = pref * hall_norm(c[k])
            trunc = pref * truncated_norm(c[k], k, Nn)
            r1 = k * math.exp(-2 * t * k * (Nn - k))
            mc_rows.append((Nn, t, k, mean, se, full, trunc, r1))
print("     N      t   k    MC mean      SE        formula     z    | r=1 only   z  | full Hall (k>N)  z")
for Nn, t, k, mean, se, full, trunc, r1 in mc_rows:
    form = trunc   # = full for k <= N
    extra = f"  {full:.5f} {(mean - full) / se:7.1f}" if k > Nn else ""
    print(f"    {Nn:2d} {t:6.2f} {k:3d}  {mean:9.5f}  {se:.5f}   {form:9.5f} {(mean - form) / se:6.2f}   |"
          f" {r1:8.5f} {(mean - r1) / se:6.1f}{extra}")
ok = all(abs(mean - trunc) < 3 * se for Nn, t, k, mean, se, full, trunc, r1 in mc_rows)
check(f"MC agrees with Claim 2 (truncated form for k > N) within 3 SE in all {len(mc_rows)} cells", ok,
      "max |z| = %.2f" % max(abs(mean - trunc) / se for Nn, t, k, mean, se, full, trunc, r1 in mc_rows))
ok = all(abs(mean - full) > 10 * se for Nn, t, k, mean, se, full, trunc, r1 in mc_rows if k > Nn)
check("MC rejects the full Hall norm for k > N (|z| > 10)", ok)
ok = all(abs(mean - r1) > 5 * se for Nn, t, k, mean, se, full, trunc, r1 in mc_rows if k <= Nn and k >= 4 and t == 0.03)
check("MC rejects the r = 1 term alone for 4 <= k <= N, t = 0.03 (|z| > 5): the cluster terms are seen", ok)

print(f"\n(time {time.time() - T0:.0f} s)")
print(f"{NCHK} checks, {NFAIL} failed")
