"""Checks for al_theorem.md: murmurations for Atkin-Lehner operators W_Q and Atkin-Lehner eigenspaces at one squarefree level.

Run: (ulimit -v 1200000; timeout 600 nice /home/bas/praime/.venv/bin/python -I -B al_check.py) > al_check_out.txt
Sections: A class numbers and local factors; B local means (exact enumeration); C trace formula of T_p W_Q on the new
space against PARI; D densities; E discovery numerics (declared as such); F predictions files for the blind test.
Ends with 'N checks, 0 failed'.
"""
import sys, os, math, json, time, hashlib, itertools
from fractions import Fraction as Fr
import numpy as np
import cypari2
from sympy import isprime, primerange, factorint, divisors

HERE = os.path.dirname(os.path.abspath(__file__))
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
NCHK = [0, 0]


def check(name, ok, info=''):
    NCHK[0] += 1
    if not ok:
        NCHK[1] += 1
    print('%-4s %s %s' % (name, 'ok  ' if ok else 'FAIL', info), flush=True)


def H(n):  # Hurwitz class number H(n), n >= 0 (H(0) = -1/12), exact
    return Fr(str(pari.qfbhclassno(n)))


def kron(a, b):
    return int(pari.kronecker(a, b))


def omega(n):
    return len(factorint(n))


def mu(n):
    f = factorint(n)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)


def sigma1(n):
    return sum(divisors(n))


def pk(k, t2, ell):  # p_k(t, ell) = ell^(k/2-1) U_{k-2}(t/(2 sqrt ell)) as polynomial in t^2 (k even)
    # recursion P_{j+1} = t P_j - ell P_{j-1}; even/odd split: keep P_j as dict power->coef
    P = [{}, {0: 1}]  # P_1 = 0, P_2 = 1
    for j in range(2, k):
        a, b = P[-1], P[-2]
        nw = {}
        for e, c in a.items():
            nw[e + 1] = nw.get(e + 1, 0) + c
        for e, c in b.items():
            nw[e] = nw.get(e, 0) - ell * c
        P.append(nw)
    return sum(c * t2 ** (e // 2) for e, c in P[-1].items())


# ---------- local data at a prime ell of a negative discriminant Delta ----------
def local_e_chi(Delta, ell):
    """e = v_ell(conductor), chi = (Delta_0/ell) for a negative discriminant Delta."""
    v = 0; u = Delta
    while u % ell == 0:
        u //= ell; v += 1
    if ell != 2:
        if v % 2 == 0:
            return v // 2, kron(u, ell)
        return v // 2, 0
    if v % 2 == 0:
        if u % 4 == 1:
            return v // 2, kron(u, 2)
        return v // 2 - 1, 0
    return (v - 3) // 2, 0


def F_loc(e, chi, ell):
    x = Fr(1, ell)
    return sum(x ** j for j in range(e)) + x ** e / (1 - chi * x)


def G_loc(Delta, ell):  # xi_Delta(ell) * F_ell(Delta) = ell^(-e)(chi-1)/(1-chi/ell)
    e, chi = local_e_chi(Delta, ell)
    return Fr(chi - 1) / (1 - Fr(chi, ell)) / ell ** e


def xi_cls(Delta, ell):  # Martin's xi_Delta(ell) from class numbers (Delta < 0)
    if Delta % (ell * ell) == 0:
        return ell * H(-Delta // (ell * ell)) / H(-Delta) - 1
    return kron(Delta, ell) - 1


def e_gen(s, ell):  # generic local mean (ell not dividing N): note, Lemma 4
    if s % ell == 0:
        return Fr(ell * ell, ell * ell - 1)
    return 1 - Fr(1, (ell - 1) ** 2 * (ell + 1))


def e_new(s, ell):  # local mean of G_ell at ell | M = N/Q (this note, Lemma B)
    if s % ell == 0:
        return Fr(-ell, ell + 1)
    return Fr(-ell * ell, ell * ell - 1)


# ---------- Martin's trace formula (arXiv 2409.02338, Prop. 5.? 'prop:Tl-large') ----------
def trace_TW(N, Q, k, ell):
    M = N // Q
    tot = Fr(0)
    smax = int(math.isqrt(4 * ell // Q)) + 1
    for s in range(-smax, smax + 1):
        if s * s * Q > 4 * ell:
            continue
        Delta = s * s * Q * Q - 4 * Q * ell
        xi = Fr(1)
        for l in factorint(M):
            xi *= xi_cls(Delta, l)
        tot += pk(k, s * s * Q, ell) * xi * H(-Delta)
    return -tot / 2 + (sigma1(ell) * mu(M) if k == 2 else 0) - (1 if N == 1 else 0)


def sectionA():
    print('== A. class numbers and local factors')
    # A1 own count of H(D) (reduced forms, weights 1/2, 1/3) against PARI
    def H_own(D):
        if D % 4 in (1, 2):
            return Fr(0)
        tot = Fr(0); b = D % 2
        while b * b <= D // 3:
            ac = (b * b + D) // 4
            a = max(b, 1)
            while a * a <= ac:
                if ac % a == 0:
                    c = ac // a
                    if a == b or a == c or b == 0:
                        w = Fr(1)
                        if a == b == c:
                            w = Fr(1, 3)
                        elif b == 0 and a == c:
                            w = Fr(1, 2)
                        tot += w
                    else:
                        tot += 2
                a += 1
            b += 2
        return tot
    bad = [D for D in range(3, 1501) if H_own(D) != H(D)]
    check('A1', not bad, 'own count of reduced forms = qfbhclassno for 3 <= D <= 1500 (%d bad)' % len(bad))
    # A2 xi_Delta(ell) F_ell(Delta) = ell^(-e)(chi-1)/(1-chi/ell) with xi from class numbers (Martin's definition)
    cnt = 0; bad = 0
    for D in range(3, 4001):
        if D % 4 in (1, 2):
            continue
        for l in (2, 3, 5, 7):
            e, chi = local_e_chi(-D, l)
            lhs = xi_cls(-D, l) * F_loc(e, chi, l)
            cnt += 1; bad += lhs != G_loc(-D, l)
            # Martin's explicit formula for ell^2 | Delta
            if D % (l * l) == 0:
                mart = Fr((l - 1) * (chi - 1), (l ** (e + 1) - 1) - chi * (l ** e - 1))
                bad += mart != xi_cls(-D, l)
    check('A2', bad == 0, '%d cases (D <= 4000, ell = 2,3,5,7): xi*F = ell^-e (chi-1)/(1-chi/ell); Martin explicit xi = class-number xi' % cnt)
    # A3 primes of Q divide the fundamental discriminant exactly (Lemma 3 with Q in place of N)
    cnt = 0; bad = 0
    for N in [6, 10, 14, 15, 21, 30, 42, 66, 70, 105, 210, 2 * 4999, 3 * 3329, 1155, 2310]:
        for Q in divisors(N):
            if Q == 1:
                continue
            for p in primerange(2, 300):
                if N % p == 0:
                    continue
                s = 0
                while s * s * Q < 4 * p:
                    D = Q * (4 * p - s * s * Q)
                    d0 = int(pari.coredisc(-D)); f0 = math.isqrt(D // abs(d0))
                    for l in factorint(Q):
                        cnt += 1; bad += not (d0 % l == 0 and f0 % l != 0)
                    s += 1
    check('A3', bad == 0, '%d cases: every prime of Q divides Delta_0 and not f_0 (D = Q(4p - s^2 Q))' % cnt)


def sectionB():
    print('== B. local means at ell | M = N/Q (exact enumeration with enclosure)')
    J = {2: 16, 3: 10, 5: 7, 7: 6}
    cnt = 0; bad = 0; wmax = 0.0
    for l in (2, 3, 5, 7):
        mod = l ** J[l]
        ps = np.array([p for p in range(1, mod) if p % l], dtype=np.int64)
        for Q in [q for q in (1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 101, 5003, 3331, 10009) if q % l][:9]:
            for s in range(0, 9):
                Dl = (Q * (s * s * Q - 4 * ps)) % mod  # Delta mod ell^J (Delta = Q(s^2 Q - 4p) < 0)
                tot = Fr(0); und = 0
                vals, counts = np.unique(Dl, return_counts=True)
                for d, c in zip(vals.tolist(), counts.tolist()):
                    if d == 0:
                        und += c; continue
                    v = 0; u = d
                    while u % l == 0:
                        u //= l; v += 1
                    need = 1 if l != 2 else 3
                    if J[l] - v < need:
                        und += c; continue
                    # representative negative discriminant with the same ell-adic data: d - mod*K
                    rep = d - mod * (10 ** 6)  # negative, same ell-adic data as Delta up to ell^J
                    tot += c * G_loc(rep, l)
                n = len(ps)
                lo = tot / n - Fr(2 * und, n); hi = tot / n + Fr(2 * und, n)
                pred = e_new(s, l)
                cnt += 1; bad += not (lo <= pred <= hi); wmax = max(wmax, float(hi - lo))
    check('B1', bad == 0, '%d cases (ell = 2,3,5,7; Q units at ell; s = 0..8): mean of xi*F over p in Z_ell^x = -ell^2/(ell^2-1) (ell !| s), -ell/(ell+1) (ell | s); enclosure width <= %.1e' % (cnt, wmax))


def sectionC():
    print('== C. trace of T_p W_Q on S_k^new(N): Martin\'s class-number formula against PARI')
    cnt = 0; bad = 0
    for N in (6, 10, 14, 15, 21, 30, 42):
        for k in (2, 4):
            tf = pari('mfcoefs(mftraceform([%d,%d],0),%d)' % (N, k, N * 150))
            for Q in divisors(N):
                for p in primerange(2, 151):
                    if N % p == 0:
                        continue
                    lhs = trace_TW(N, Q, k, p)
                    rhs = (-1) ** omega(Q) * Fr(Q) ** (1 - k // 2) * Fr(str(tf[Q * p]))
                    cnt += 1; bad += lhs != rhs
    check('C1', bad == 0, '%d triples (N in 6,10,14,15,21,30,42; k = 2,4; all Q | N; primes p <= 150): exact' % cnt)
    # C2 trace-form route against the Atkin-Lehner matrices of PARI (independent of a_f(q) = -q^(k/2-1) w_q)
    cnt = 0; bad = 0
    for N in (30, 42, 70):
        for k in (2, 4):
            pari('mf=mfinit([%d,%d],0)' % (N, k))
            for Q in divisors(N):
                if Q == 1:
                    continue
                pari('A=mfatkininit(mf,%d); W=A[2]*A[3]' % Q)
                for p in primerange(2, 60):
                    if N % p == 0:
                        continue
                    rhs = Fr(str(pari('trace(W*mfheckemat(mf,%d))' % p)))
                    cnt += 1; bad += trace_TW(N, Q, k, p) != rhs
    check('C2', bad == 0, '%d triples (N = 30, 42, 70; k = 2, 4; Q > 1; p < 60): formula = trace(W_Q T_p) from mfatkininit' % cnt)
    # C3 traces of W_Q (ell = 1) and the eigenspace decomposition against PARI orbits
    cnt = 0; bad = 0
    pari('trc(c,d) = if(type(c)=="t_POLMOD", trace(c), d*c)')
    for N in (30, 42, 66, 70, 78, 105):
        for k in (2, 4):
            pari('mf=mfinit([%d,%d],0); L=mfeigenbasis(mf); F=mffields(mf)' % (N, k))
            norb = int(pari('#L'))
            if norb == 0:
                continue
            primesN = sorted(factorint(N))
            signs = {}
            for q in primesN:
                signs[q] = [int(x[0]) for x in pari('mfatkineigenvalues(mf,%d)' % q)]
            degs = [int(pari('poldegree(F[%d])' % (i + 1))) for i in range(norb)]
            for Q in divisors(N):
                if Q == 1:
                    continue
                direct = sum(degs[i] * math.prod(signs[q][i] for q in factorint(Q)) for i in range(norb))
                cnt += 1; bad += trace_TW(N, Q, k, 1) != direct
            for p in (5, 7, 11, 13):
                if N % p == 0:
                    continue
                ap = [Fr(str(pari('trc(mfcoef(L[%d],%d),%d)' % (i + 1, p, degs[i])))) for i in range(norb)]
                for eps in itertools.product((1, -1), repeat=len(primesN)):
                    ed = dict(zip(primesN, eps))
                    direct = sum(ap[i] for i in range(norb) if all(signs[q][i] == ed[q] for q in primesN))
                    comb = sum(math.prod(ed[q] for q in factorint(Q)) * (trace_TW(N, Q, k, p)) for Q in divisors(N)) / 2 ** len(primesN)
                    cnt += 1; bad += comb != direct
    check('C3', bad == 0, '%d cases (N in 30,42,66,70,78,105; k = 2,4): tr W_Q (ell = 1) and tr(T_p | eigenspace) = 2^-w sum_Q eps(Q) tr(T_p W_Q) against PARI orbits' % cnt)
    # C4 exact splitting for N = m q, Q0 | m: t_Q0 at level mq = Leg - t_Q0 at level m  (q^2 !| Delta)
    cnt = 0; bad = 0
    for m, q in ((2, 41), (2, 53), (3, 37), (6, 31), (5, 43), (2, 101)):
        N = m * q
        for k in (2, 4):
            for Q0 in divisors(m):
                for p in primerange(2, 400):
                    if N % p == 0 or 4 * p * Q0 >= q * q:
                        continue
                    leg = Fr(0); smax = int(math.isqrt(4 * p // Q0)) + 1
                    for s in range(-smax, smax + 1):
                        if s * s * Q0 > 4 * p:
                            continue
                        Delta = s * s * Q0 * Q0 - 4 * Q0 * p
                        xi = Fr(1)
                        for l in factorint(m // Q0):
                            xi *= xi_cls(Delta, l)
                        leg += pk(k, s * s * Q0, p) * xi * kron(Delta, q) * H(-Delta)
                    leg = -leg / 2
                    lhs = trace_TW(N, Q0, k, p)
                    rhs = leg - trace_TW(m, Q0, k, p) if m > 1 or Q0 > 1 else None
                    cnt += 1; bad += lhs != rhs
    check('C4', bad == 0, '%d cases: tr(T_p W_Q0 | S^new(mq)) = Leg_Q0(p) - tr(T_p W_Q0 | S^new(m)) exactly (Q0 | m, 4 p Q0 < q^2)' % cnt)


# ---------- densities ----------
C0 = None


def C0_val():
    global C0
    if C0 is None:
        ls = np.array(list(primerange(2, 10 ** 6)), dtype=float)
        C0 = float(np.exp(np.sum(np.log1p(-1.0 / ((ls - 1) ** 2 * (ls + 1))))))
    return C0


def nu_prime(s):
    v = C0_val()
    for l in factorint(s):
        v *= float(e_gen(0, l) / e_gen(1, l))
    return v


def Evals(N, Q, smax):
    """E^(Q,M)(s), s = 0..smax: mean of the local factors (Lemma 3', Lemma 4, Lemma B)."""
    M = N // Q
    out = []
    for s in range(smax + 1):
        if s == 0:
            v = math.pi ** 2 / 6
            for l in factorint(Q):
                v *= 1 - 1 / l ** 2
            for l in factorint(M):
                v *= -(l - 1) / l
        else:
            v = nu_prime(s)
            for l in factorint(N):
                v /= float(e_gen(s, l))
            for l in factorint(M):
                v *= float(e_new(s, l))
        out.append(v)
    return out


def U(kk, x):
    if kk == 0:
        return 1.0
    if kk == 2:
        return 4 * x * x - 1
    if kk == 4:
        return 16 * x ** 4 - 12 * x * x + 1
    raise ValueError


def That(N, Q, k, p, E):
    """model of t_Q(p) = p^(1-k/2) tr(T_p W_Q | S_k^new(N)): class numbers replaced by sqrt(D)/pi * local mean."""
    M = N // Q
    tot = 0.0; s = 0
    while s * s * Q < 4 * p:
        D = 4 * p * Q - s * s * Q * Q
        term = U(k - 2, s * math.sqrt(Q) / (2 * math.sqrt(p))) * math.sqrt(D) / math.pi * E[s]
        tot += term if s == 0 else 2 * term
        s += 1
    return -tot / 2 + (mu(M) * (p + 1) if k == 2 else 0)


def density(N, Q, k, y):
    """D_k^(Q,M)(y) at y = p/N: the limit of t_Q(p)/dim with dim -> (k-1) phi(N)/12."""
    p = y * N
    E = Evals(N, Q, int(math.isqrt(int(4 * p / Q) + 1)) + 1)
    phi = N
    for l in factorint(N):
        phi = phi // l * (l - 1)
    return That(N, Q, k, p, E) / ((k - 1) * phi / 12) - (mu(N // Q) * 12 / ((k - 1) * phi) if k == 2 else 0)


def M_note(N, k, y):  # formula (2.1) of murmurations_one_level.md
    kap = N; pi_ = 1.0
    for l in factorint(N):
        kap = kap * l / (l - 1); pi_ *= 1 + 1 / l
    kap /= N
    tot = 2 * math.pi * pi_ * math.sqrt(y); r = 1
    while r * r < 4 * y:
        nu = nu_prime(r)
        for l in factorint(N):
            nu /= float(e_gen(r, l))
        tot += (-1) ** (k // 2 - 1) * 12 / math.pi * kap * nu * U(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    return tot / (k - 1) - (12 * kap * y if k == 2 else 0)


def sectionD():
    print('== D. densities')
    # D1 the W_N density is (-1)^(k/2) M^(N) of the note
    mx = 0.0
    for N in (2 * 10007, 3 * 6673, 30, 1000003):
        for k in (2, 4):
            for y in np.linspace(0.01, 4, 200):
                mx = max(mx, abs(density(N, N, k, y) - (-1) ** (k // 2) * M_note(N, k, y)))
    check('D1', mx < 1e-9, 'D^(N,1) = (-1)^(k/2) M^(N) (formula (2.1)) on 1600 points, max diff %.1e' % mx)
    # D2 s = 0 term: (-1)^(k/2+w(M)) (2 pi/(k-1)) pi(Q) sqrt(y/M) for y M < 1/4
    mx = 0.0
    for N, Q in ((2 * 10007, 10007), (3 * 6673, 6673), (30, 5), (30, 15), (42, 7)):
        M = N // Q
        for k in (2, 4):
            for y in np.linspace(0.001, 0.249 / M, 50):
                piQ = math.prod(1 + 1 / l for l in factorint(Q))
                kap = N / math.prod(l - 1 for l in factorint(N)) / math.prod(1 for _ in [0])
                kap = N
                for l in factorint(N):
                    kap = kap * (l - 1) / l
                kap = N / kap
                pred = (-1) ** (k // 2 + omega(M)) * 2 * math.pi / (k - 1) * piQ * math.sqrt(y / M) + (12 * kap * mu(M) * y if k == 2 else 0)
                mx = max(mx, abs(density(N, Q, k, y) - pred))
    check('D2', mx < 1e-9, 'below y = 1/(4M) the W_Q density is (-1)^(k/2+w(M)) (2pi/(k-1)) pi(Q) sqrt(y/M) + 12 kappa(N) mu(M) y [k=2]; max diff %.1e' % mx)
    # D3 table, q -> infinity (q = 10^9+7)
    q = 1000000007
    print('     y      ' + '  '.join('%9s' % s for s in ('D_q(2q)', 'D_N(2q)', '++', '+-', '-+', '--', 'D_q(3q)', 'D_N(3q)', '++', '+-', '-+', '--')))
    for y in (0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0):
        row = []
        for m in (2, 3):
            dq = density(m * q, q, 2, y); dn = density(m * q, m * q, 2, y)
            row += [dq, dn] + [eq * dq + e1 * eq * dn for e1 in (1, -1) for eq in (1, -1)]
        print('    %4.2f   ' % y + '  '.join('%9.4f' % v for v in row))
    check('D3', True, 'table printed (k = 2; eigenspace (e1, eq) density = eq D_q + e1 eq D_N; columns ++,+-,-+,-- = (e_p1, e_q))')
    # D4 Martin's Figure (S_2(2q), 3000 < q < 6000): level-averaged eigenspace densities M^eps(x; beta = 2), X = 6000
    print('     Martin Fig. pq-profile3000 (levels 6000..12000), predicted M^eps(l/X; 2) and M^eps/4:')
    for ell in (250, 800, 2000, 3700, 4000, 5500):
        vals = []
        ts = np.linspace(1, 2, 401)
        for e1, eq in ((1, 1), (-1, -1), (1, -1), (-1, 1)):
            f = [(eq * density(2 * q, q, 2, ell / (6000 * t)) + e1 * eq * density(2 * q, 2 * q, 2, ell / (6000 * t))) * t for t in ts]
            vals.append(np.trapezoid(f, ts) / 1.5)
        print('     l = %4d: (++,--,+-,-+) = %s ; /4 = %s' % (ell, ' '.join('%6.3f' % v for v in vals), ' '.join('%6.3f' % (v / 4) for v in vals)))
    check('D4', True, 'Martin figure comparison printed (readings of the figure: see note)')




# ---------- E. discovery numerics (weight 2, levels N = p1 q) ----------
WIN = [(0.25 * i, 0.25 * (i + 1)) for i in range(16)]
EPS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]  # (sign at p1, sign at q)


def htable(nmax):
    return np.array([-1.0 / 12] + [float(x) for x in pari('vector(%d,n,qfbhclassno(n))' % nmax)])


def kron_small(D, l):  # Kronecker (D/l) for numpy int arrays, l = 2 or odd prime (small)
    if l == 2:
        r = np.mod(D, 8)
        return np.where(r % 2 == 0, 0, np.where((r == 1) | (r == 7), 1, -1))
    r = np.mod(D, l)
    qr = np.zeros(l, dtype=np.int64)
    for x in range(1, l):
        qr[x * x % l] = 1
    return np.where(r == 0, 0, np.where(qr[r] == 1, 1, -1))


def xiH_np(D, l, Ht):  # xi_D(l) H(|D|) for negative discriminants D (numpy), via the table
    A = -D
    div = (A % (l * l) == 0)
    out = (kron_small(D, l) - 1) * Ht[A]
    if div.any():
        out = np.where(div, l * Ht[np.where(div, A // (l * l), 0)] - Ht[A], out)
    return out


def level_traces(N, p1, q, Ht, hcl):
    """t_Q(p) for Q = 1, p1, q, N at every prime p <= 4N, p !| N (weight 2; Martin's formula)."""
    ps = np.array([p for p in primerange(2, 4 * N + 1) if N % p], dtype=np.int64)
    leg = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        leg[x * x % q] = 1
    legq = lambda D: np.where(np.mod(D, q) == 0, 0, np.where(leg[np.mod(D, q)] == 1, 1, -1))
    out = {}
    for Q0 in (1, p1):  # small discriminants: table
        smax = np.floor(np.sqrt(4 * ps / Q0)).astype(np.int64) + 1
        idx = np.repeat(np.arange(len(ps)), smax + 1)
        s = np.concatenate([np.arange(m + 1) for m in smax])
        P = ps[idx]
        ok = s * s * Q0 < 4 * P
        idx, s, P = idx[ok], s[ok], P[ok]
        D = s * s * Q0 * Q0 - 4 * Q0 * P
        val = (legq(D) - 1) * (xiH_np(D, p1, Ht) if Q0 == 1 else Ht[-D])
        w = np.where(s == 0, 0.5, 1.0)
        tot = np.bincount(idx, weights=w * val, minlength=len(ps))
        out[Q0] = -tot + (ps + 1) * (1 if Q0 == 1 else -1)   # mu(N) = +1, mu(q) = -1
    for Q in (q, N):  # large discriminants: PARI
        M = N // Q
        res = np.zeros(len(ps))
        for i, p in enumerate(ps.tolist()):
            tot = 0.0; s = 0
            while s * s * Q < 4 * p:
                A = 4 * p * Q - s * s * Q * Q
                if M == 1:
                    v = float(hcl(A))
                else:
                    if A % (M * M) == 0:
                        v = M * float(hcl(A // (M * M))) - float(hcl(A))
                    else:
                        v = (kron(-A, M) - 1) * float(hcl(A))
                tot += 0.5 * v if s == 0 else v
                s += 1
            res[i] = -tot + mu(M) * (p + 1)
        out[Q] = res
    return ps, out


def level_models(N, p1, q, ps):
    EQ = Evals(N, q, 8); EN = Evals(N, N, 8)
    TQ = np.array([That(N, q, 2, p, EQ) for p in ps.tolist()])
    TN = np.array([That(N, N, 2, p, EN) for p in ps.tolist()])
    return TQ, TN


def eps_dims(N, p1, q):
    dim = int(pari('mfdim([%d,2],0)' % N))
    tw = {Q: trace_TW(N, Q, 2, 1) for Q in (p1, q, N)}
    de = [(dim + e1 * tw[p1] + eq * tw[q] + e1 * eq * tw[N]) / 4 for e1, eq in EPS]
    return dim, [int(x) for x in de], {str(k): str(v) for k, v in tw.items()}


def wmeans(ps, N, X):
    """log-weighted window means and single-value standard errors of X (rows: quantities)."""
    y = ps / N; w = np.log(ps)
    mean = np.zeros((X.shape[0], 16)); se = np.zeros((X.shape[0], 16)); cnt = []
    for i, (a, b) in enumerate(WIN):
        m = (y > a) & (y <= b); ww = w[m]; cnt.append(int(m.sum()))
        for j in range(X.shape[0]):
            x = X[j][m]; mu_ = np.sum(ww * x) / np.sum(ww)
            mean[j, i] = mu_; se[j, i] = math.sqrt(np.sum(ww ** 2 * (x - mu_) ** 2)) / np.sum(ww)
    return mean, se, cnt


DISC = {'2q': [(2, 5003), (2, 5009), (2, 5011)], '3q': [(3, 3331), (3, 3343), (3, 3347)]}


def sectionE():
    print('== E. discovery numerics (declared: seen data), weight 2, all p <= 4N, 16 windows of y = p/N')
    hcl = pari('D -> qfbhclassno(D)')
    # E0: vectorised data against the exact formula (validated in C1-C3) at small and discovery levels
    mx = 0.0; cnt = 0
    for p1, q in ((2, 101), (3, 67), (2, 5003)):
        N = p1 * q
        Ht = htable(16 * p1 * N + 8)
        ps, t = level_traces(N, p1, q, Ht, hcl)
        sel = [i for i, p in enumerate(ps.tolist()) if p < 60 or i % max(1, len(ps) // 25) == 0]
        for i in sel:
            for Q in (1, p1, q, N):
                mx = max(mx, abs(t[Q][i] - float(trace_TW(N, Q, 2, int(ps[i]))))); cnt += 1
    check('E0', mx < 1e-6, '%d values t_Q(p) (N = 202, 201, 10006; all Q): vectorised = exact formula, max diff %.1e' % (cnt, mx))
    disc = {}
    for T, lv in DISC.items():
        for p1, q in lv:
            N = p1 * q; t0 = time.time()
            Ht = htable(16 * p1 * N + 8)
            ps, t = level_traces(N, p1, q, Ht, hcl)
            TQ, TN = level_models(N, p1, q, ps)
            dim, de, tw = eps_dims(N, p1, q)
            rows = []
            for k_, (e1, eq) in enumerate(EPS):
                rows.append((t[1] + e1 * t[p1] + eq * t[q] + e1 * eq * t[N]) / (4 * de[k_]))  # data E^eps(p)
            for k_, (e1, eq) in enumerate(EPS):
                rows.append((eq * TQ + e1 * eq * TN) / (4 * de[k_]))  # own model
            for k_, (e1, eq) in enumerate(EPS):
                rows.append(e1 * eq * TN / (4 * de[k_]))  # W_N-only model
            for e1 in (1, -1):
                rows.append((t[1] + e1 * t[p1]) / dim)  # symmetric parts (expected -> 0, Lemma T)
                rows.append((t[q] + e1 * t[N]) / dim)  # antisymmetric parts (unconditional theorem)
                rows.append((TQ + e1 * TN) / dim)
            for e1 in (1, -1):
                rows.append(e1 * TN / dim)
            mean, se, cnt = wmeans(ps, N, np.array(rows))
            disc[N] = {'type': T, 'p1': p1, 'q': q, 'dim': dim, 'dimeps': de, 'trW': tw, 'nprimes': cnt,
                       'data': mean[0:4].tolist(), 'se': se[0:4].tolist(), 'own': mean[4:8].tolist(), 'wn': mean[8:12].tolist(),
                       'sym': [mean[12].tolist(), mean[15].tolist()], 'anti': [mean[13].tolist(), mean[16].tolist()],
                       'anti_model': [mean[14].tolist(), mean[17].tolist()], 'sym_se': [se[12].tolist(), se[15].tolist()],
                       'anti_se': [se[13].tolist(), se[16].tolist()], 'anti_wn': [mean[18].tolist(), mean[19].tolist()]}
            print('     N = %d = %d*%d: dim %d, eigenspace dims %s, %d primes, %.0f s' % (N, p1, q, dim, de, len(ps), time.time() - t0), flush=True)
    json.dump(disc, open(os.path.join(HERE, 'discovery_data.json'), 'w'))
    stats_and_small(disc, hcl)


def chi2_block(disc, Ns, sig, sig_a):
    """chi2 of the eigenspace means (own model and alternatives) and of the antisymmetric combinations."""
    out = {'own': 0.0, 'wn': 0.0, 'eps0': 0.0, 'eps1': 0.0, 'eps2': 0.0, 'eps3': 0.0, 'anti_own': 0.0, 'anti_wn': 0.0}
    for N in Ns:
        d = disc[N]
        for k_ in range(4):
            x = np.array(d['data'][k_])
            out['own'] += np.sum(((x - np.array(d['own'][k_])) / sig[k_]) ** 2)
            out['wn'] += np.sum(((x - np.array(d['wn'][k_])) / sig[k_]) ** 2)
            for j in range(4):
                if j != k_:
                    alt = np.array(d['own'][j]) * d['dimeps'][j] / d['dimeps'][k_]
                    out['eps%d' % j] += np.sum(((x - alt) / sig[k_]) ** 2)
        for i in range(2):
            x = np.array(d['anti'][i])
            out['anti_own'] += np.sum(((x - np.array(d['anti_model'][i])) / sig_a[i]) ** 2)
            out['anti_wn'] += np.sum(((x - np.array(d['anti_wn'][i])) / sig_a[i]) ** 2)
    out['alt_min'] = min(out[k_] for k_ in ('wn', 'eps0', 'eps1', 'eps2', 'eps3') if out[k_] > 0)
    return out


def errors_from(disc, Ns):
    """frozen errors: per type one scale c_T (between-level scatter of the residuals / single-value SE) x the rms SE."""
    R = np.array([[np.array(disc[N]['data'][k_]) - np.array(disc[N]['own'][k_]) for k_ in range(4)] for N in Ns])
    SE = np.sqrt(np.mean(np.array([disc[N]['se'] for N in Ns]) ** 2, axis=0))
    cT = math.sqrt(np.sum((R - R.mean(axis=0)) ** 2) / (len(Ns) - 1) / np.sum(SE ** 2))
    Ra = np.array([[np.array(disc[N]['anti'][i]) - np.array(disc[N]['anti_model'][i]) for i in range(2)] for N in Ns])
    SEa = np.sqrt(np.mean(np.array([disc[N]['anti_se'] for N in Ns]) ** 2, axis=0))
    ca = math.sqrt(np.sum((Ra - Ra.mean(axis=0)) ** 2) / (len(Ns) - 1) / np.sum(SEa ** 2))
    return cT * SE, ca * SEa, cT, ca


def stats_and_small(disc, hcl):
    summary = {}
    for T in DISC:
        Ns = [p1 * q for p1, q in DISC[T]]
        sig, sig_a, cT, ca = errors_from(disc, Ns)
        c2 = chi2_block(disc, Ns, sig, sig_a)
        S = np.array([[disc[N]['sym'][i] for N in Ns] for i in range(2)])  # e1 x L x w
        corr = [float(np.mean([np.corrcoef(S[i][a], S[i][b])[0, 1] for a in range(3) for b in range(a + 1, 3)])) for i in range(2)]
        summary[T] = {'sig': sig.tolist(), 'sig_anti': sig_a.tolist(), 'cT': cT, 'c_anti': ca, 'chi2': c2,
                      'sym_rms': float(np.sqrt(np.mean(S ** 2))), 'sym_corr': corr}
        print('     type %s (levels %s): errors c_T = %.2f x SE (eigenspaces), %.2f x SE (antisymmetric parts)' % (T, Ns, cT, ca))
        print('       eigenspaces, 192 terms: chi2 own %.1f; W_N-only %.1f; other eigenspace densities %s' % (
            c2['own'], c2['wn'], ', '.join('%s %.0f' % (EPS[j], c2['eps%d' % j]) for j in range(4))))
        print('       antisymmetric parts (t_q + e1 t_N)/dim, 96 terms: chi2 own %.1f; W_N-only %.1f' % (c2['anti_own'], c2['anti_wn']))
        print('       symmetric parts (t_1 + e1 t_p1)/dim: rms of window means %.3f; mean correlation between levels (e1 = +1, -1): %.2f, %.2f' % (
            summary[T]['sym_rms'], corr[0], corr[1]))
        for k_ in range(4):
            print('       eps %-8s data : %s' % (EPS[k_], ' '.join('%5.2f' % np.mean([disc[N]['data'][k_][w] for N in Ns]) for w in range(16))))
            print('       eps %-8s model: %s' % (EPS[k_], ' '.join('%5.2f' % np.mean([disc[N]['own'][k_][w] for N in Ns]) for w in range(16))))
    json.dump(summary, open(os.path.join(HERE, 'discovery_summary.json'), 'w'))
    check('E1', all(summary[T]['chi2']['own'] <= 384 and summary[T]['chi2']['anti_own'] <= 192 for T in summary),
          'in-sample: chi2 own <= 2 per term (eigenspaces 192 terms, antisymmetric parts 96 terms) for both types')
    check('E2', all(summary[T]['chi2']['alt_min'] >= summary[T]['chi2']['own'] + 300 and summary[T]['chi2']['anti_wn'] >= summary[T]['chi2']['anti_own'] + 150 for T in summary),
          'every alternative worse by >= 300 (eigenspaces) and the W_N-only antisymmetric model by >= 150')
    # E3 the symmetric parts at smaller levels (decay with N; evidence for Lemma T, not a proof)
    for T, lv in (('2q', [(2, 503), (2, 509), (2, 521)]), ('3q', [(3, 337), (3, 347), (3, 349)])):
        rms = []
        for p1, q in lv:
            N = p1 * q
            Ht = htable(16 * p1 * N + 8)
            ps, t = level_traces(N, p1, q, Ht, hcl)
            dim = int(pari('mfdim([%d,2],0)' % N))
            X = np.array([(t[1] + e1 * t[p1]) / dim for e1 in (1, -1)])
            mean, se, cnt = wmeans(ps, N, X)
            rms.append(float(np.sqrt(np.mean(mean ** 2))))
        print('     symmetric parts at N ~ 1e3 (%s): rms %.3f (N ~ 1e4: %.3f; ratio %.2f, N^(-1/2) scaling would give %.2f)' % (
            T, float(np.sqrt(np.mean(np.array(rms) ** 2))), summary[T]['sym_rms'], summary[T]['sym_rms'] / float(np.sqrt(np.mean(np.array(rms) ** 2))), math.sqrt(1000 / 10000)))
    check('E3', True, 'decay of the symmetric parts printed (reported, not a check of a claim)')




# ---------- F. predictions for the blind test (no data at any test level is computed here) ----------
def blind_levels():
    def nxt(lo, n, excl):
        out = []; q = int(lo) + 1
        while len(out) < n:
            if isprime(q) and q not in excl:
                out.append(q)
            q += 1
        return out
    lv = {'2q': [(2, q) for q in nxt(20000 / 2, 10, {10007})], '3q': [(3, q) for q in nxt(20000 / 3, 10, {6673})]}
    forbidden = set(range(30030, 30531)) | {2 * 10007, 3 * 6673, 6 * 3343, 5 * 4003, 9998, 9987, 9978}
    assert all(p1 * q not in forbidden for T in lv for p1, q in lv[T])
    return lv


def model_windows(p1, q):
    N = p1 * q
    ps = np.array([p for p in primerange(2, 4 * N + 1) if N % p], dtype=np.int64)
    TQ, TN = level_models(N, p1, q, ps)
    dim, de, tw = eps_dims(N, p1, q)
    rows = [(eq * TQ + e1 * eq * TN) / (4 * de[k_]) for k_, (e1, eq) in enumerate(EPS)]
    rows += [e1 * eq * TN / (4 * de[k_]) for k_, (e1, eq) in enumerate(EPS)]
    rows += [(TQ + e1 * TN) / dim for e1 in (1, -1)] + [e1 * TN / dim for e1 in (1, -1)]
    mean, se, cnt = wmeans(ps, N, np.array(rows))
    own = mean[0:4].tolist()
    alt = [[(mean[j] * de[j] / de[k_]).tolist() for j in range(4)] for k_ in range(4)]  # alt[k][j]: density of eigenspace j, normalised as k
    return {'p1': p1, 'q': q, 'dim': dim, 'dimeps': de, 'trW': tw, 'nprimes': cnt, 'own': own, 'wn': mean[4:8].tolist(), 'alt': alt,
            'anti_model': mean[8:10].tolist(), 'anti_wn': mean[10:12].tolist()}


def sectionF():
    print('== F. predictions (frozen before any test-level data)')
    summ = json.load(open(os.path.join(HERE, 'discovery_summary.json')))
    dec = {'D1_own_chi2_per_term_max': 2.0, 'D2_alt_delta_per_level_min': 100.0, 'D2_alt_chi2_per_term_min': 2.0,
           'D3_anti_own_chi2_per_term_max': 2.0, 'D3_anti_wn_delta_per_level_min': 50.0}
    for which, lv in (('blind', blind_levels()), ('dry', DISC)):
        doc = {'meta': {'windows': WIN, 'eps_order': EPS, 'weight': 2, 'decision': dec,
                        'sigma': {T: summ[T]['sig'] for T in DISC}, 'sigma_anti': {T: summ[T]['sig_anti'] for T in DISC},
                        'statement': 'eigenspace (e_p1, e_q) mean of a_f(p) over primes p in windows of y = p/N -> e_q D^(q,p1) + e_p1 e_q D^(N,1) '
                                     '(al_theorem.md, Theorem B; the symmetric part (t_1 + e_p1 t_p1)/dim -> 0 is Lemma T, under GRH); '
                                     'antisymmetric parts (t_q + e1 t_N)/dim -> D^(q,p1) + e1 D^(N,1) (Theorem A, no hypothesis)'},
               'levels': {T: [p1 * q for p1, q in lv[T]] for T in lv}, 'pred': {}}
        for T in lv:
            for p1, q in lv[T]:
                doc['pred'][str(p1 * q)] = model_windows(p1, q)
        txt = json.dumps(doc, sort_keys=True).encode()
        fn = os.path.join(HERE, 'predictions_%s.json' % which)
        open(fn, 'wb').write(txt)
        print('     %s: %s  sha256 %s' % (which, {T: doc['levels'][T] for T in doc['levels']}, hashlib.sha256(txt).hexdigest()))
        # expected separations from the predictions alone (no data)
        for T in lv:
            sig = np.array(summ[T]['sig']); siga = np.array(summ[T]['sig_anti'])
            sep = {'wn': 0.0, 'eps': 1e18, 'anti_wn': 0.0}
            epsj = [0.0] * 4
            for p1, q in lv[T]:
                P = doc['pred'][str(p1 * q)]
                for k_ in range(4):
                    sep['wn'] += np.sum(((np.array(P['own'][k_]) - np.array(P['wn'][k_])) / sig[k_]) ** 2)
                    for j in range(4):
                        if j != k_:
                            epsj[j] += np.sum(((np.array(P['own'][k_]) - np.array(P['alt'][k_][j])) / sig[k_]) ** 2)
                for i in range(2):
                    sep['anti_wn'] += np.sum(((np.array(P['anti_model'][i]) - np.array(P['anti_wn'][i])) / siga[i]) ** 2)
            print('       type %s: expected separation (sum (own - alt)^2 / sigma^2): W_N-only %.0f; other eigenspaces %s; antisymmetric W_N-only %.0f' % (
                T, sep['wn'], ', '.join('%.0f' % x for x in epsj), sep['anti_wn']))
    check('F1', True, 'predictions written (blind levels: none of 30030-30530, 2*10007, 3*6673, 6*3343, 5*4003, 9998, 9987, 9978)')
    import re
    reg = dict(re.findall(r"'(blind|dry)': '([0-9a-f]{64})'", open(os.path.join(HERE, 'al_blind_test.py')).read()))
    got = {w: hashlib.sha256(open(os.path.join(HERE, 'predictions_%s.json' % w), 'rb').read()).hexdigest() for w in ('blind', 'dry')}
    check('F2', got == reg, 'regenerated predictions files have the sha256 registered in al_blind_test.py (byte-for-byte reproducible)')


if __name__ == '__main__':
    t0 = time.time()
    parts = sys.argv[1:] or ['A', 'B', 'C', 'D', 'E', 'F']
    for prt in parts:
        globals()['section' + prt]()
    print('%d checks, %d failed  (%.0f s)' % (NCHK[0], NCHK[1], time.time() - t0))
