#!/usr/bin/env python3
"""martin_check.py -- checks for martin_theorem.md (Martin, arXiv 2409.02338, Conjecture 1.10, case r = m+1).

Own code. Run: (ulimit -v 1200000; timeout 600 nice python -I -B martin_check.py)
Sections
  A  trace formula (Martin Prop. 5.5) against PARI trace forms; eigenspace dimensions against PARI;
     the splitting t_{Q0}(Pq) = Leg - t_{Q0}(P); t_{Q0}(2) = 0 at weight 2 (S_2(2) = 0)
  B  local Fubini: the q-mean at fixed l equals the l-mean at fixed q (exact histograms mod v^J);
     local means of F_v and G_v over q at fixed l against e_v(s), e_v^new(s) (rigorous enclosures)
  C  class-number means over primes q at fixed l against E(s) (key lemma (a), empirical)
  D  Lemma T in Martin's order: multiplicity of fundamental discriminants among Delta_s <= 6 tau(4 l Q0)
  F  densities D_q, D_N for 2q against the table of al_theorem.md
  E  A^eps(l; 10^4, 2*10^4) for N = 2q, k = 2, twelve primes l, against the finite-X model and M^eps(x; 2)
Numerics use PARI's qfbhclassno; for arguments above 5*10^5 PARI's class-group routine may assume GRH
(computational caveat only; no proof step uses these numbers).
"""
import math, sys, time
from fractions import Fraction as Fr
import numpy as np
import gmpy2
import cypari2

pari = cypari2.Pari()
pari.allocatemem(400 * 10**6)
T0 = time.time()
NCHK = 0
NFAIL = 0


def check(name, ok, info=""):
    global NCHK, NFAIL
    NCHK += 1
    if not ok:
        NFAIL += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {info}", flush=True)


_H = {}


def H(n):
    """Hurwitz class number (exact Fraction); H(0) not used."""
    n = int(n)
    if n % 4 in (1, 2):
        return Fr(0)
    if n not in _H:
        _H[n] = Fr(str(pari.qfbhclassno(n)))
    return _H[n]


def Hf(n):
    return float(H(n))


def kron(a, n):
    return int(gmpy2.kronecker(int(a), int(n)))


def primefac(n):
    return [int(p) for p in pari.factor(n)[0]]


def mu(n):
    f = pari.factor(n)
    return 0 if any(int(e) > 1 for e in f[1]) else (-1) ** len(f[0])


def xi_v(Delta, v):
    """Martin's xi_Delta(v) for a prime v, Delta < 0 a discriminant (exact)."""
    if Delta % (v * v) != 0:
        return Fr(kron(Delta, v) - 1)
    return v * H(-Delta // (v * v)) / H(-Delta) - 1


def xi_M(Delta, M):
    r = Fr(1)
    for v in (primefac(M) if M > 1 else []):
        r *= xi_v(Delta, v)
    return r


def U(n, x):
    """Chebyshev U_n(x) for exact rational x^2 handled by caller; here numeric or Fraction via recursion."""
    a, b = 1, 2 * x
    if n == 0:
        return a
    for _ in range(n - 1):
        a, b = b, 2 * x * b - a
    return b


def U_even_sq(n, x2):
    """U_n(x) for even n as a polynomial in x^2 (exact for rational x^2)."""
    # U_0 = 1, U_2 = 4x^2 - 1, U_4 = 16x^4 - 12x^2 + 1, U_6 = 64x^6 - 80x^4 + 24x^2 - 1
    tab = {0: [1], 2: [-1, 4], 4: [1, -12, 16], 6: [-1, 24, -80, 64]}
    c = tab[n]
    return sum(ci * x2 ** i for i, ci in enumerate(c))


def t_Q(l, N, Q, k):
    """l^(1-k/2) tr(T_l W_Q | S_k^new(N)) by Martin Prop. 5.5 (exact Fraction), l prime, N squarefree."""
    M = N // Q
    tot = Fr(0)
    s = 0
    while s * s * Q < 4 * l:
        Delta = s * s * Q * Q - 4 * Q * l
        w = 1 if s == 0 else 2
        tot += w * U_even_sq(k - 2, Fr(s * s * Q, 4 * l)) * xi_M(Delta, M) * H(-Delta)
        s += 1
    val = -tot / 2
    if N == 1:
        val -= Fr(1) if k == 2 else Fr(1, l ** (k // 2 - 1))
    if k == 2:
        val += mu(M) * (l + 1)
    return val


def trW_Q(N, Q, k):
    """tr(W_Q | S_k^new(N)) for Q > 1 by Prop. 5.5 with l = 1 (exact)."""
    M = N // Q
    tot = Fr(0)
    s = 0
    while s * s * Q <= 4:
        Delta = s * s * Q * Q - 4 * Q
        w = 1 if s == 0 else 2
        tot += w * U_even_sq(k - 2, Fr(s * s * Q, 4)) * xi_M(Delta, M) * H(-Delta)
        s += 1
    val = -tot / 2
    if k == 2:
        val += mu(M)
    return val


def primes_in(a, b):
    return [int(p) for p in pari.primes([a, b])]


# ---------------------------------------------------------------- A: trace formula and eigenspaces
print("== A: trace formula (Prop. 5.5) against PARI ==", flush=True)
levels = [10, 14, 22, 26, 34, 38, 46, 58, 15, 21, 33, 39, 51, 57, 30, 42, 66, 78]
nA1 = 0
badA1 = []
for N in levels:
    for k in (2, 4):
        if int(pari(f"mfdim([{N},{k}],0)")) == 0:
            continue
        tf = pari(f"mftraceform([{N},{k}],0)")
        divs = [int(d) for d in pari.divisors(N)]
        cmax = max(divs) * 47
        co = [Fr(str(c)) for c in pari.mfcoefs(tf, cmax)]
        for l in primes_in(2, 47):
            if N % l == 0:
                continue
            for Q in divs:
                om = len(primefac(Q)) if Q > 1 else 0
                tr_pari = (-1) ** om * Fr(1, Q ** (k // 2 - 1)) * co[Q * l]
                tr_mart = t_Q(l, N, Q, k) * Fr(l ** (k // 2 - 1))
                nA1 += 1
                if tr_pari != tr_mart:
                    badA1.append((N, k, l, Q, tr_pari, tr_mart))
check("A1 Prop. 5.5 = PARI trace forms, tr(T_l W_Q) exact", not badA1,
      f"({nA1} triples N,k,l,Q; 18 levels 2q,3q,6q; k=2,4; l<=47; failures {badA1[:2]})")

# eigenspace dimensions and eigenspace traces against PARI newforms with AL eigenvalues
nA2 = 0
badA2 = []
for N in [26, 38, 46, 58, 21, 33, 39, 57, 42, 66, 78]:
    for k in (2, 4):
        if int(pari(f"mfdim([{N},{k}],0)")) == 0:
            continue
        mf = pari(f"mfinit([{N},{k}],0)")
        ps = primefac(N)
        alv = [pari.mfatkineigenvalues(mf, v) for v in ps]
        norb = len(alv[0])
        B = pari.mfeigenbasis(mf)
        dims = {}
        traces = {}
        for o in range(norb):
            sig = tuple(int(alv[i][o][0]) for i in range(len(ps)))
            deg = len(alv[0][o])
            dims[sig] = dims.get(sig, 0) + deg
            co = pari.mfcoefs(B[o], 50)
            # trace over the orbit of a_l: trace of the coefficient field element
            for l in primes_in(2, 47):
                if N % l == 0:
                    continue
                c = co[l]
                trc = (Fr(str(pari.trace(c))) if str(pari.type(c)) == 't_POLMOD' else deg * Fr(str(c)))
                traces[(sig, l)] = traces.get((sig, l), Fr(0)) + trc
        import itertools
        divs = [int(d) for d in pari.divisors(N)]
        dimnew = int(pari(f"mfdim([{N},{k}],0)"))
        for sig in itertools.product((1, -1), repeat=len(ps)):
            eps = dict(zip(ps, sig))
            epsQ = lambda Q: math.prod(eps[v] for v in primefac(Q)) if Q > 1 else 1
            d_form = (dimnew + sum(epsQ(Q) * trW_Q(N, Q, k) for Q in divs if Q > 1)) / 2 ** len(ps)
            nA2 += 1
            if d_form != dims.get(sig, 0):
                badA2.append((N, k, sig, d_form, dims.get(sig, 0)))
            for l in primes_in(2, 47):
                if N % l == 0:
                    continue
                tr_form = sum(epsQ(Q) * t_Q(l, N, Q, k) for Q in divs) / 2 ** len(ps) * Fr(l ** (k // 2 - 1))
                nA2 += 1
                if tr_form != traces.get((sig, l), Fr(0)):
                    badA2.append((N, k, sig, l, tr_form, traces.get((sig, l))))
check("A2 eigenspace dims and eigenspace traces 2^-r sum_Q eps(Q) t_Q = PARI newforms", not badA2,
      f"({nA2} cases; failures {badA2[:2]})")

# A3: splitting t_{Q0}(Pq) = Leg - t_{Q0}(P) (exact), and A4: t_{Q0}(P) = 0 for P = 2, 3, k = 2
nA3 = 0
badA3 = []
for P in (2, 3, 6):
    for q in primes_in(11, 60):
        if P % q == 0:
            continue
        N = P * q
        for Q0 in [int(d) for d in pari.divisors(P)]:
            for k in (2, 4):
                for l in primes_in(2, 120):
                    if N % l == 0 or q * q <= 4 * l * Q0:
                        continue
                    M0 = P // Q0
                    leg = Fr(0)
                    s = 0
                    while s * s * Q0 < 4 * l:
                        Delta = s * s * Q0 * Q0 - 4 * Q0 * l
                        w = 1 if s == 0 else 2
                        leg += w * U_even_sq(k - 2, Fr(s * s * Q0, 4 * l)) * xi_M(Delta, M0) * kron(Delta, q) * H(-Delta)
                        s += 1
                    leg = -leg / 2
                    nA3 += 1
                    if t_Q(l, N, Q0, k) != leg - t_Q(l, P, Q0, k):
                        badA3.append((P, q, Q0, k, l))
check("A3 t_Q0(Pq) = Leg_Q0 - t_Q0(P) exactly (q^2 > 4 l Q0)", not badA3, f"({nA3} cases; failures {badA3[:3]})")
nA4 = 0
badA4 = []
for P in (2, 3):
    for Q0 in [int(d) for d in pari.divisors(P)]:
        for l in primes_in(2, 3000):
            if P % l == 0:
                continue
            nA4 += 1
            if t_Q(l, P, Q0, 2) != 0:
                badA4.append((P, Q0, l))
check("A4 t_Q0(P) = 0 for P = 2, 3, k = 2 (S_2(P) = 0: Hurwitz-type class number relation)", not badA4,
      f"({nA4} cases)")

# ---------------------------------------------------------------- B: local Fubini (exact)
print("== B: local means over q at fixed l ==", flush=True)


def canon_table(v, J):
    mod = v ** J
    sq = sorted({(c * c) % mod for c in range(mod) if c % v})
    r = np.arange(mod, dtype=np.int64)
    best = r.copy()
    for c2 in sq:
        best = np.minimum(best, (r * c2) % mod)
    return best


okB1 = True
nB1 = 0
for v, J in ((2, 7), (3, 5), (5, 3), (7, 3)):
    mod = v ** J
    can = canon_table(v, J)
    units = np.array([u for u in range(1, mod) if u % v], dtype=np.int64)
    for Q0 in (1, 3, 5, 7, 11, 15):
        if Q0 % v == 0:
            continue
        for s in range(0, 7):
            # D(l, q) = Q0 q (4 l - s^2 Q0 q)
            Dmat = (Q0 * units[None, :] % mod) * ((4 * units[:, None] - s * s * Q0 * units[None, :]) % mod) % mod
            C = can[Dmat]  # rows: l, cols: q
            ref = np.bincount(C[0], minlength=mod)
            for i in range(len(units)):
                nB1 += 2
                if not np.array_equal(np.bincount(C[i], minlength=mod), ref):
                    okB1 = False
                if not np.array_equal(np.bincount(C[:, i], minlength=mod), ref):
                    okB1 = False
check("B1 homogeneity: distribution of D(l,q) mod v^J up to unit squares, over q, is the same for every unit l "
      "and equals the distribution over l at fixed q", okB1, f"({nB1} histogram comparisons; v = 2,3,5,7)")


def local_FG(v, D, J):
    """(F_v, G_v) of D > 0 known mod v^J (as an integer residue); None if undetermined."""
    mod = v ** J
    if v != 2:
        d = D % mod
        if d == 0:
            return None
        e = 0
        while d % v == 0:
            d //= v
            e += 1
        if e % 2 == 1:
            chi, a = 0, (e - 1) // 2
        else:
            chi, a = kron(-d, v), e // 2
    else:
        w = (-D) % mod
        bits = J
        a = 0
        while True:
            if bits < 3:
                return None
            if w % 4 == 1:
                chi = 1 if w % 8 == 1 else -1
                break
            if w % 4 in (2, 3):
                raise ValueError("not a discriminant")
            if bits < 4:
                return None
            if (w // 4) % 4 in (2, 3):
                chi = 0
                break
            w //= 4
            bits -= 2
            a += 1
    x = 1.0 / v
    F = sum(x ** j for j in range(a)) + x ** a / (1 - chi * x)
    G = x ** a * (chi - 1) / (1 - chi * x)
    return F, G


def e_loc(v, s):
    return v * v / (v * v - 1) if s % v == 0 else 1 - 1 / ((v - 1) ** 2 * (v + 1))


def e_new(v, s):
    return -v * v / (v * v - 1) if s % v else -v / (v + 1)


okB2 = True
nB2 = 0
maxw = 0
for v, J in ((2, 13), (3, 8), (5, 6), (7, 5)):
    mod = v ** J
    units = [u for u in range(1, mod) if u % v]
    for l in [1, 3, 5, 7][: (4 if v == 2 else 2)] + [mod - 1]:
        if l % v == 0:
            continue
        for Q0 in (1, 3, 5, 7):
            if Q0 % v == 0:
                continue
            for s in range(0, 5):
                sF = sG = 0.0
                und = 0
                for q in units:
                    D = Q0 * q * (4 * l - s * s * Q0 * q)
                    r = local_FG(v, D, J)
                    if r is None:
                        und += 1
                    else:
                        sF += r[0]
                        sG += r[1]
                n = len(units)
                Flo, Fhi = sF / n, (sF + und * v / (v - 1)) / n
                Glo, Ghi = (sG - 2 * und) / n, sG / n
                maxw = max(maxw, Fhi - Flo, Ghi - Glo)
                nB2 += 2
                if not (Flo - 1e-12 <= e_loc(v, s) <= Fhi + 1e-12):
                    okB2 = False
                    print("  F miss", v, l, Q0, s, Flo, Fhi, e_loc(v, s))
                if not (Glo - 1e-12 <= e_new(v, s) <= Ghi + 1e-12):
                    okB2 = False
                    print("  G miss", v, l, Q0, s, Glo, Ghi, e_new(v, s))
check("B2 q-means at fixed l: mean F_v = e_v(s), mean G_v = e_v^new(s) inside rigorous enclosures", okB2,
      f"({nB2} cases, max width {maxw:.1e})")
# v | Q0: F_v = 1 for all q (v odd and v = 2)
okB3 = True
nB3 = 0
for v, Q0 in ((3, 3), (5, 5), (3, 15), (2, 2), (2, 6), (7, 21)):
    for l in range(1, 200, 2):
        if l % v == 0 or math.gcd(l, Q0) > 1:
            continue
        for q in range(1, 200, 2):
            if q % v == 0 or math.gcd(q, Q0) > 1:
                continue
            for s in range(0, 5):
                D = Q0 * q * (4 * l - s * s * Q0 * q)
                r = local_FG(v, D, 12)
                nB3 += 1
                if r is None or abs(r[0] - 1) > 1e-12:
                    okB3 = False
check("B3 local factor F_v = 1 at v | Q0 (v exactly divides the fundamental discriminant), every q", okB3, f"({nB3} cases)")

# B4: expansion of q -> chi_v(q (a - c q)) in characters psi mod v: coefficients have modulus sqrt(v)/(v-1)
# (psi not in {psi0, chi_v}) or 1/(v-1); for v | c it is chi_v(a) * chi_v(q) exactly
okB4 = True
nB4 = 0
for v in (3, 5, 7, 11, 13, 17):
    g = int(pari.znprimroot(v).lift())
    logt = {pow(g, i, v): i for i in range(v - 1)}
    for a in range(1, v):
        for c in range(0, v):
            vals = [kron(x * (a - c * x), v) for x in range(1, v)]
            for j in range(v - 1):  # psi_j(g^i) = exp(2 pi i j i/(v-1))
                coef = sum(vals[x - 1] * complex(math.cos(-2 * math.pi * j * logt[x] / (v - 1)),
                                                 math.sin(-2 * math.pi * j * logt[x] / (v - 1))) for x in range(1, v)) / (v - 1)
                nB4 += 1
                if c % v == 0:
                    want = 1.0 if j == (v - 1) // 2 else 0.0
                elif j in (0, (v - 1) // 2):
                    want = 1 / (v - 1)
                else:
                    want = math.sqrt(v) / (v - 1)
                if abs(abs(coef) - want) > 1e-9:
                    okB4 = False
check("B4 Jacobi-sum expansion of chi_v(q(4l - s^2 Q0 q)): |a(psi)| = sqrt(v)/(v-1) or 1/(v-1); pure character if v | s",
      okB4, f"({nB4} coefficients, v <= 17)")

# ---------------------------------------------------------------- densities
C0 = 1.0
for v in primes_in(2, 10**6):
    C0 *= 1 - 1 / ((v - 1) ** 2 * (v + 1))


from functools import lru_cache


@lru_cache(maxsize=None)
def nuprime(s):
    r = C0
    for v in primefac(s):
        r *= e_loc(v, s) / (1 - 1 / ((v - 1) ** 2 * (v + 1)))
    return r


ZETA2 = math.pi ** 2 / 6


@lru_cache(maxsize=None)
def _Efac(s, Qp, Mp, qbig):
    Qp, Mp = list(Qp), list(Mp)
    """E^(Q,M)(s): Qp = primes of Q (excluding the big prime), Mp = primes of M; qbig: the big prime of Q or None (limit)."""
    if s == 0:
        r = ZETA2
        for v in Qp:
            r *= 1 - v ** -2
        if qbig:
            r *= 1 - qbig ** -2
        for v in Mp:
            r *= -(v - 1) / v
        return r
    r = nuprime(s)
    for v in Qp + Mp:
        r /= e_loc(v, s)
    if qbig:
        r /= e_loc(qbig, s)
    for v in Mp:
        r *= e_new(v, s)
    return r


def Efac(s, Qp, Mp, qbig=None):
    return _Efac(s, tuple(Qp), tuple(Mp), qbig)


def kappa(ps):
    r = 1.0
    for v in ps:
        r *= v / (v - 1)
    return r


def Dens(y, Qp, Mp, k=2):
    """q -> infinity limit of D_k^(Q0 q, M)(y), (1.1) of al_theorem.md, y = l/N."""
    M = math.prod(Mp) if Mp else 1
    kap = kappa(Qp + Mp)
    val = (-1) ** (k // 2) * Efac(0, Qp, Mp) * math.sqrt(y / M)
    s = 1
    while s * s < 4 * y * M:
        x2 = s * s / (4 * y * M)
        val -= (1 / M) * U_even_sq(k - 2, x2) * Efac(s, Qp, Mp) * math.sqrt(4 * y * M - s * s)
        s += 1
    val *= 12 * kap / ((k - 1) * math.pi)
    if k == 2:
        val += 12 * kap * (-1) ** len(Mp) * y
    return val


# ---------------------------------------------------------------- F: densities against al_theorem table
print("== F: densities ==", flush=True)
ys = [0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 3]
tab_Dq = [-0.207, -0.995, 0.827, -0.719, -2.425, 1.440, -0.458, -2.774]
tab_DN = [-0.907, -0.580, 0.585, -1.713, 2.366, -1.273, 1.924, -2.291]
tab_3q = [-0.089, -0.022, -0.310, 0.705, -0.094, 0.612, 0.116, -1.012]
dev = max(max(abs(Dens(y, [], [2]) - a), abs(Dens(y, [2], []) - b), abs(Dens(y, [], [3]) - c))
          for y, a, b, c in zip(ys, tab_Dq, tab_DN, tab_3q))
check("F1 D_q (2q), D_N (2q), D_q (3q) = table of al_theorem.md (3 decimals)", dev < 1.5e-3, f"(max dev {dev:.1e})")


def Mlim(x, Qp_list, signs, beta=2.0, n=4000):
    """M(x; beta) = int_1^beta sum_i sign_i D_i(x/t) t dt / int_1^beta t dt (midpoint rule)."""
    if x == 0:
        return 0.0
    h = (beta - 1) / n
    tot = 0.0
    for i in range(n):
        t = 1 + (i + 0.5) * h
        tot += sum(sg * Dens(x / t, Qp, Mp) for sg, (Qp, Mp) in zip(signs, Qp_list)) * t * h
    return tot / ((beta ** 2 - 1) / 2)


# continuity at 0 and in x (spot): M(x) -> 0 as x -> 0
m_small = [abs(Mlim(x, [([], [2]), ([2], [])], [1, 1], n=400)) for x in (1e-4, 1e-3)]
check("F2 M^eps(x; 2) -> 0 as x -> 0 (eps = ++)", m_small[0] < 0.05 and m_small[0] < m_small[1], f"({m_small})")

# ---------------------------------------------------------------- C: class-number means over primes q
print("== C: class-number means over primes q at fixed l (key lemma (a)) ==", flush=True)
l0 = int(pari.nextprime(20000))
zs = []
for P, qa, qb in ((2, 5000, 10000), (3, 3334, 6666)):
    qs = [q for q in primes_in(qa, qb) if q != l0]
    for Qp, Mp in (([], [P]), ([P], [])):
        Q0 = math.prod(Qp) if Qp else 1
        M0 = math.prod(Mp) if Mp else 1
        for s in (0, 1, 2):
            vals = []
            for q in qs:
                Q = Q0 * q
                if s * s * Q >= 4 * l0:
                    continue
                D = Q * (4 * l0 - s * s * Q)
                Delta = -D
                vals.append(float(xi_M(Delta, M0)) * math.pi * Hf(D) / math.sqrt(D))
            vals = np.array(vals)
            pred = np.mean([Efac(s, Qp, Mp, qbig=q) for q in qs[:5]])
            z = (vals.mean() - pred) / (vals.std(ddof=1) / math.sqrt(len(vals)))
            zs.append(z)
            print(f"   P={P} Q={'Q0' if Qp else ''}q s={s}: n={len(vals)} mean {vals.mean():+.4f} pred E(s) {pred:+.4f} z {z:+.2f}")
check("C1 mean over primes q of xi_Delta(M) L*(D_s(q)) = E(s) at l = 20011 (|z| <= 4 each)", max(abs(z) for z in zs) <= 4,
      f"(12 cases, max |z| {max(abs(z) for z in zs):.2f})")

# ---------------------------------------------------------------- D: multiplicities of fundamental discriminants
print("== D: Lemma T in Martin's order ==", flush=True)
X = 10**4
xs = [0.15, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5, 1.8, 2.1, 2.5, 2.9]
ells = [int(pari.nextprime(int(x * X))) for x in xs]
okD = True
maxfib = 0
for l in ells:
    for Q0 in (1, 2, 3, 6):
        fib = {}
        s = 0
        while s * s * Q0 < 4 * l:
            Delta = s * s * Q0 * Q0 - 4 * Q0 * l
            d0 = int(pari.coredisc(Delta))
            fib[d0] = fib.get(d0, 0) + 1
            s += 1
        bound = 6 * int(pari.numdiv(4 * l * Q0))
        maxfib = max(maxfib, max(fib.values()))
        if max(fib.values()) > bound:
            okD = False
check("D1 each fundamental discriminant occurs for at most 6 tau(4 l Q0) values s >= 0", okD,
      f"(12 primes l, Q0 = 1,2,3,6; max fibre {maxfib})")

# ---------------------------------------------------------------- E: Martin's average for N = 2q
print("== E: A^eps(l; 10^4, 2*10^4), N = 2q, k = 2 ==", flush=True)
# Decision rule fixed before the run: E1 chi2(own model, 48 terms, errors pooled per l over the four eigenspaces
# from three interleaved groups of levels) <= 100; E2 every alternative (W_N only, the three other eigenspace models,
# zero) has chi2 >= chi2_own + 100; E3 rms of the Lemma-T part <= 0.25 * rms of the model.
qs = primes_in(5000, 10000)
assert 2 * qs[0] >= X and 2 * qs[-1] <= 2 * X
dimnew = {q: int(pari(f"mfdim([{2*q},2],0)")) for q in qs}
trW2 = {q: -Fr(1, 2) * ((kron(-8, q) - 1) * H(8) + 2 * (kron(-4, q) - 1) * H(4)) - 1 for q in qs}
trWq = {q: -Fr(1, 2) * xi_v(-4 * q, 2) * H(4 * q) - 1 for q in qs}
trW2q = {q: -Fr(1, 2) * H(8 * q) + 1 for q in qs}
EPS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]  # (eps_2, eps_q)
groups = [qs[i::3] for i in range(3)]
rows = []
for x, l in zip(xs, ells):
    # q-independent pieces for Q = 1, 2
    S1 = []
    s = 0
    while s * s < 4 * l:
        Delta = s * s - 4 * l
        S1.append((Delta, (1 if s == 0 else 2) * float(xi_v(Delta, 2)) * Hf(-Delta)))
        s += 1
    S2 = []
    s = 0
    while 2 * s * s < 4 * l:
        Delta = 4 * s * s - 8 * l
        S2.append((Delta, (1 if s == 0 else 2) * Hf(-Delta)))
        s += 1
    per = {}
    for q in qs:
        den = {e: float((dimnew[q] + e[0] * trW2[q] + e[1] * trWq[q] + e[0] * e[1] * trW2q[q]) / 4) for e in EPS}
        if q == l:
            per[q] = (den, None, None, None)
            continue
        t1 = -0.5 * sum(w * (kron(Dl, q) - 1) for Dl, w in S1) + (l + 1)
        t2 = -0.5 * sum(w * (kron(Dl, q) - 1) for Dl, w in S2) - (l + 1)
        tq = 0.0
        mq = 0.0
        s = 0
        while s * s * q < 4 * l:
            Delta = s * s * q * q - 4 * q * l
            D = -Delta
            w = 1 if s == 0 else 2
            tq += -0.5 * w * float(xi_v(Delta, 2)) * Hf(D)
            mq += -0.5 * w * math.sqrt(D) / math.pi * Efac(s, [], [2], qbig=q)
            s += 1
        tq -= (l + 1)
        mq -= (l + 1)
        t2q = 0.0
        m2q = 0.0
        s = 0
        while 2 * s * s * q < 4 * l:
            D = 8 * q * l - 4 * s * s * q * q
            w = 1 if s == 0 else 2
            t2q += -0.5 * w * Hf(D)
            m2q += -0.5 * w * math.sqrt(D) / math.pi * Efac(s, [2], [], qbig=q)
            s += 1
        t2q += (l + 1)
        m2q += (l + 1)
        per[q] = (den, (t1, t2, tq, t2q), (mq, m2q), None)
    for e in EPS:
        e2, eq = e
        def agg(qlist):
            num = sum((per[q][1][0] + e2 * per[q][1][1] + eq * per[q][1][2] + e2 * eq * per[q][1][3]) / 4
                      for q in qlist if per[q][1] is not None)
            den = sum(per[q][0][e] for q in qlist)
            mq_ = sum(eq * per[q][2][0] / 4 for q in qlist if per[q][1] is not None)
            mN_ = sum(e2 * eq * per[q][2][1] / 4 for q in qlist if per[q][1] is not None)
            lt = sum((per[q][1][0] + e2 * per[q][1][1]) / 4 for q in qlist if per[q][1] is not None)
            return num / den, mq_ / den, mN_ / den, lt / den
        A, mqp, mNp, LT = agg(qs)
        Ag = [agg(g)[0] for g in groups]
        se = float(np.std(Ag, ddof=1) / math.sqrt(3))
        Mx = eq * Mlim(l / X, [([], [2]), ([2], [])], [1, e2], n=600)
        rows.append(dict(x=l / X, l=l, e=e, A=A, model=mqp + mNp, modelN=mNp, M=Mx, LT=LT, se=se))
    print(f"   l={l} done ({time.time()-T0:.0f}s)", flush=True)

# pooled errors per l over the four eigenspaces
for l in ells:
    rr = [r for r in rows if r["l"] == l]
    sp = math.sqrt(np.mean([r["se"] ** 2 for r in rr]))
    for r in rr:
        r["sp"] = sp
print("   x       l     eps     A^eps    model_X   M^eps(x;2)  LemmaT-part   se(pooled)")
for r in rows:
    print(f"  {r['x']:.4f} {r['l']:6d} {r['e'][0]:+d}{r['e'][1]:+d}  {r['A']:+8.3f}  {r['model']:+8.3f}  {r['M']:+8.3f}"
          f"     {r['LT']:+7.3f}     {r['sp']:.3f}")


def chi2(key):
    return sum((r["A"] - key(r)) ** 2 / r["sp"] ** 2 for r in rows)


byle = {(r["l"], r["e"]): r for r in rows}
c_own = chi2(lambda r: r["model"])
calt = {"W_N only": chi2(lambda r: r["modelN"]), "zero": chi2(lambda r: 0.0)}
for sh in EPS[1:]:
    # the model of the eigenspace (e2*sh2, eq*sh1): every other eigenspace, at the same l
    calt[f"other eigenspace {sh}"] = chi2(lambda r, sh=sh: byle[(r["l"], (r["e"][0] * sh[0], r["e"][1] * sh[1]))]["model"])
print(f"   chi2 (48 terms): own model {c_own:.1f}; " + "; ".join(f"{k} {v:.1f}" for k, v in calt.items()))
c_lim = chi2(lambda r: r["M"])
print(f"   chi2 against the limit M^eps(x; 2): {c_lim:.1f}; max |model_X - M^eps| "
      f"{max(abs(r['model'] - r['M']) for r in rows):.3f}")
check("E1 chi2(own finite-X model) <= 100 on 48 terms", c_own <= 100, f"({c_own:.1f})")
worst = min(calt.values())
check("E2 every alternative (W_N only, the three other eigenspace models, zero) has chi2 >= own + 100",
      worst >= c_own + 100, f"(smallest alternative {worst:.1f})")
rms_lt = math.sqrt(np.mean([r["LT"] ** 2 for r in rows]))
rms_mod = math.sqrt(np.mean([r["model"] ** 2 for r in rows]))
check("E3 rms of the Lemma-T part (Q = 1, 2 terms) <= 0.25 rms(model)", rms_lt <= 0.25 * rms_mod,
      f"(rms {rms_lt:.3f} vs model rms {rms_mod:.3f})")
print(f"   rms(A - model_X) {math.sqrt(np.mean([(r['A']-r['model'])**2 for r in rows])):.3f}; "
      f"rms(A - M^eps) {math.sqrt(np.mean([(r['A']-r['M'])**2 for r in rows])):.3f}; "
      f"rms pooled se {math.sqrt(np.mean([r['sp']**2 for r in rows])):.3f}")

bias = {e: sum(float((dimnew[q] + e[0] * trW2[q] + e[1] * trWq[q] + e[0] * e[1] * trW2q[q]) / 4) for q in qs)
        / (sum(dimnew[q] for q in qs) / 4) for e in EPS}
print("   eigenspace size ratio sum|H^eps| / (sum dim / 4): " + ", ".join(f"{e}: {b:.3f}" for e, b in bias.items()))


def lemmaT_rms(Xv):
    qv = primes_in(Xv // 2, Xv)
    dn = {q: int(pari(f"mfdim([{2*q},2],0)")) for q in qv}
    out = []
    for x in xs:
        l = int(pari.nextprime(int(x * Xv)))
        S1 = []
        s = 0
        while s * s < 4 * l:
            Delta = s * s - 4 * l
            S1.append((Delta, (1 if s == 0 else 2) * float(xi_v(Delta, 2)) * Hf(-Delta)))
            s += 1
        S2 = []
        s = 0
        while 2 * s * s < 4 * l:
            Delta = 4 * s * s - 8 * l
            S2.append((Delta, (1 if s == 0 else 2) * Hf(-Delta)))
            s += 1
        t1s = t2s = 0.0
        for q in qv:
            if q == l:
                continue
            t1s += -0.5 * sum(w * (kron(Dl, q) - 1) for Dl, w in S1) + (l + 1)
            t2s += -0.5 * sum(w * (kron(Dl, q) - 1) for Dl, w in S2) - (l + 1)
        den = sum(dn.values()) / 4
        for e2 in (1, -1):
            out.append((t1s + e2 * t2s) / 4 / den)
    return math.sqrt(np.mean(np.square(out)))


lt_small = lemmaT_rms(2500)
lt_big = lemmaT_rms(10000)
print(f"   Lemma-T part (Q = 1, 2), rms over the 12 x and both signs at 2: X = 2500: {lt_small:.3f}; X = 10^4: {lt_big:.3f};"
      f" ratio {lt_big/lt_small:.2f} (X^-1/4: 0.71, X^-1/2: 0.50, X^-3/4: 0.35)")

# ---------------------------------------------------------------- G: Martin's figure (illustration only)
print("== G: M^eps(x; 2)/4 at the points of the figure S_2(2q), 3000 < q < 6000 (X = 6000) ==", flush=True)
for lf in (250, 800, 2000, 4000):
    vals = [eq * Mlim(lf / 6000, [([], [2]), ([2], [])], [1, e2], n=600) / 4 for e2, eq in [(1, 1), (-1, -1), (1, -1), (-1, 1)]]
    print(f"   l={lf}: (++, --, +-, -+) M/4 = " + ", ".join(f"{v:+.2f}" for v in vals))

print(f"time {time.time()-T0:.0f}s")
print(f"{NCHK} checks, {NFAIL} failed")
