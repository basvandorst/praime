"""uncond_check.py -- numerical checks for mprime_unconditional.md (own code).

Groups (each assertion is one numbered check):
  A  identities: Hurwitz = sum over m of L(1, Kronecker) / m; periodicity of the terms T(m,n;p);
     Jacobi-sum coefficients of x -> chi_q(4x - c) on (Z/q)^x.
  B  local and global main term: sum over (m,n) of e(m,n)/(mn) = E(r) = nu'(r)(1+O(N^-2)), r = 0..6;
     the 'principal n only' sum differs for r >= 1.
  C  truncation of the Dirichlet series per discriminant (N = 10007): error vs M; mean over p of L*.
  D  ranges (ii)-(iv) at N = 10007, J = (N, 2N]: shifted-prime sums S_q = sum log p chi_q(4p - r^2 N),
     the large-sieve chain on real data, the block sums against the bound of step U5.
Run: (ulimit -v 1200000; timeout 600 nice python -I -B uncond_check.py)
"""
import math, sys, time
import numpy as np
import gmpy2
from scipy.special import psi as digamma
import cypari2

pari = cypari2.Pari()
T0 = time.time()
NCHK = 0
NFAIL = 0


def check(name, ok, info=""):
    global NCHK, NFAIL
    NCHK += 1
    if not ok:
        NFAIL += 1
    print(f"[{'ok' if ok else 'FAIL'}] {name}: {info}", flush=True)


# ---------------------------------------------------------------- tools
def sieve_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def phi(n):
    r = n
    for l in factor(n):
        r = r // l * (l - 1)
    return r


def jacobi_vec(a, n):
    """Jacobi symbol (a/n) for int64 arrays, n odd positive."""
    a = np.array(a, dtype=np.int64) % np.array(n, dtype=np.int64)
    n = np.array(n, dtype=np.int64).copy()
    a = np.broadcast_to(a, np.broadcast(a, n).shape).copy()
    n = np.broadcast_to(n, a.shape).copy()
    res = np.ones(a.shape, dtype=np.int64)
    act = a != 0
    while act.any():
        while True:
            ev = act & (a % 2 == 0)
            if not ev.any():
                break
            a[ev] //= 2
            r8 = n % 8
            res[ev & ((r8 == 3) | (r8 == 5))] *= -1
        idx = np.nonzero(act)[0]
        aa = a[idx]; nn = n[idx]
        fl = (aa % 4 == 3) & (nn % 4 == 3)
        res[idx[fl]] *= -1
        a[idx] = nn % aa
        n[idx] = aa
        act = a != 0
    return np.where(n == 1, res, 0)


def kron2(D):
    if D % 2 == 0:
        return 0
    return 1 if D % 8 in (1, 7) else -1


def kron_vec(D, n):
    """Kronecker (D/n) for a scalar D and an array of n >= 1."""
    n = np.asarray(n, dtype=np.int64)
    v = np.zeros_like(n); m = n.copy()
    while True:
        e = (m % 2 == 0)
        if not e.any():
            break
        v[e] += 1; m[e] //= 2
    k2 = kron2(D)
    out = jacobi_vec(D % m, m)
    if k2 == 0:
        out = np.where(v > 0, 0, out)
    elif k2 == -1:
        out = out * np.where(v % 2 == 1, -1, 1)
    return out


def disc_ms(D):
    """m >= 1 with m^2 | D and -D/m^2 = 0,1 mod 4 (m runs over the divisors of the square part)."""
    s = math.isqrt(D // int(pari.core(D)))
    return [int(m) for m in pari.divisors(s) if (-D // (int(m) ** 2)) % 4 in (0, 1)]


# ---------------------------------------------------------------- A: identities
print("== A: identities")
rng = np.random.default_rng(1)
aa = rng.integers(-10**6, 10**6, 4000); nn = 2 * rng.integers(0, 50000, 4000) + 1
jv = jacobi_vec(aa, nn)
check("A0 vectorised Jacobi = gmpy2.jacobi (4000 random pairs)",
      all(int(jv[i]) == gmpy2.jacobi(int(aa[i]), int(nn[i])) for i in range(4000)))
nn2 = rng.integers(1, 10**5, 3000); DD = -rng.integers(3, 10**6, 3000)
ok = True
for i in range(3000):
    D = int(DD[i])
    if D % 4 not in (0, 1):
        D = 4 * (D // 4)
    ok &= int(kron_vec(D, np.array([nn2[i]]))[0]) == gmpy2.kronecker(D, int(nn2[i]))
check("A0b vectorised Kronecker = gmpy2.kronecker (random discriminants)", ok)


def hurwitz_brute(D):
    """Hurwitz class number of D by counting reduced forms (all, not only primitive)."""
    H = 0.0
    a = 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b * b + D) % (4 * a):
                continue
            c = (b * b + D) // (4 * a)
            if c < a:
                continue
            if c == a and b < 0:
                continue
            w = 1.0
            if a == c and b == 0:
                w = 0.5
            if a == b == c:
                w = 1.0 / 3
            H += w
        a += 1
    return H


def L1_kron(Delta):
    """L(1, (Delta/.)) exactly (double precision) via the digamma formula; Delta < 0 a discriminant."""
    q = -Delta
    a = np.arange(1, q + 1)
    chi = kron_vec(Delta, a).astype(float)
    return -float(np.sum(chi * digamma(a / q))) / q


maxerr = 0.0; maxerr2 = 0.0; cnt = 0
for D in range(3, 1501):
    if D % 4 not in (0, 3):
        continue
    S = sum(L1_kron(-D // (m * m)) / m for m in disc_ms(D))
    Hp = float(pari.qfbhclassno(D))
    maxerr = max(maxerr, abs(math.sqrt(D) / math.pi * S - Hp))
    if D <= 500:
        maxerr2 = max(maxerr2, abs(hurwitz_brute(D) - Hp))
    cnt += 1
check("A1 H(D) = (sqrt D/pi) sum_m m^-1 L(1,(-D/m^2 | .)) for all D <= 1500", maxerr < 1e-8,
      f"{cnt} discriminants, max |diff| = {maxerr:.1e}")
check("A2 PARI qfbhclassno = brute-force reduced forms, D <= 500", maxerr2 < 1e-9, f"max diff {maxerr2:.1e}")


def T_term(m, n, p, r, N):
    D = N * (4 * p - r * r * N)
    if D % (m * m):
        return 0
    Dl = -D // (m * m)
    if Dl % 4 not in (0, 1):
        return 0
    return gmpy2.kronecker(Dl, n)


# periodicity of T(m, n; .) in p modulo 2 n m^2 (all integers p, sign of D arbitrary)
ok = True; tested = 0
for (N, r) in [(13, 1), (13, 2), (29, 3), (11, 0), (17, 4)]:
    for m in range(1, 7):
        for n in [1, 2, 3, 4, 5, 6, 8, 9, 12, 15, 16, 25, 27, 45]:
            P = 2 * n * m * m
            for p in range(1, 3 * P, max(1, P // 7)):
                ok &= T_term(m, n, p, r, N) == T_term(m, n, p + P, r, N) == T_term(m, n, p + 5 * P, r, N)
                tested += 1
check("A3 T(m,n;p) is periodic in p modulo 2nm^2", ok, f"{tested} tests")


# characters mod q (general q), via CRT and discrete logs
def char_table(q):
    f = factor(q)
    comps = []  # (modulus l^k, generator, order) components; residue -> exponent vector
    for l, k in f.items():
        lk = l ** k
        if l == 2:
            if k == 2:
                comps.append((lk, 3, 2))
            elif k >= 3:
                comps.append((lk, lk - 1, 2)); comps.append((lk, 5, 1 << (k - 2)))
        else:
            g = 2
            while True:
                if all(pow(g, (l - 1) // t, l) != 1 for t in factor(l - 1)) and (k == 1 or pow(g, l - 1, l * l) != 1):
                    break
                g += 1
            comps.append((lk, g, lk // l * (l - 1)))
    units = [a for a in range(q) if math.gcd(a, q) == 1]
    E = np.zeros((q, len(comps)), dtype=np.int64)
    # discrete logs per component
    bycomp = {}
    for ci, (lk, g, o) in enumerate(comps):
        bycomp.setdefault(lk, []).append(ci)
    for lk, cis in bycomp.items():
        if len(cis) == 1:
            ci = cis[0]; g, o = comps[ci][1], comps[ci][2]
            tab = {}; x = 1
            for e in range(o):
                tab[x] = e; x = x * g % lk
            for a in units:
                E[a, ci] = tab[a % lk]
        else:  # 2^k, k >= 3: a = (-1)^s 5^t
            c1, c2 = cis; o2 = comps[c2][2]
            tab = {}; x = 1
            for t in range(o2):
                tab[x] = (0, t); tab[(-x) % lk] = (1, t); x = x * 5 % lk
            for a in units:
                s, t = tab[a % lk]; E[a, c1] = s; E[a, c2] = t
    orders = np.array([c[2] for c in comps], dtype=np.int64) if comps else np.zeros(0, dtype=np.int64)
    if len(comps) == 0:
        tups = [()]
    else:
        tups = np.array(np.meshgrid(*[np.arange(o) for o in orders], indexing="ij")).reshape(len(comps), -1).T
    U = np.zeros(q, bool); U[units] = True
    chars = []
    for t in tups:
        ph = (E * (np.array(t, dtype=float) / orders if len(comps) else 0)).sum(axis=1) if len(comps) else np.zeros(q)
        v = np.exp(2j * np.pi * ph) * U
        chars.append(v)
    chars = np.array(chars)
    # primitivity
    prim = np.ones(len(chars), bool)
    for l in f:
        sub = [a for a in units if (a - 1) % (q // l) == 0]
        for i, v in enumerate(chars):
            if np.allclose(v[sub], 1.0):
                prim[i] = False
    princ = np.array([np.allclose(v[units], 1.0) for v in chars])
    return chars, prim, princ, units


def coeffs(q, c):
    """coefficients ahat(psi) of x -> chi_q(4x - c) on (Z/q)^x, chi_q = Jacobi (./q)."""
    chars, prim, princ, units = char_table(q)
    x = np.arange(q)
    f = jacobi_vec(4 * x - c, np.full(q, q)).astype(float)
    U = np.zeros(q, bool); U[units] = True
    ah = (chars.conj() * (f * U)).sum(axis=1) / len(units)
    recon = (ah[:, None] * chars).sum(axis=0)
    return ah, princ, np.max(np.abs(recon[U] - f[U]))


N0 = 10007
ok_size = True; ok_rec = True; info = []
for q in [3, 5, 7, 11, 13, 15, 21, 33, 35, 39, 55, 105, 165, 195, 231]:
    for r in [1, 2, 4, 5, 7]:
        c = r * r * N0
        if math.gcd(q, c) != 1:
            continue
        ah, princ, rec = coeffs(q, c)
        ph = phi(q)
        ok_rec &= rec < 1e-9
        a0 = abs(ah[princ][0])
        ok_size &= abs(a0 - 1 / ph) < 1e-12 and np.max(np.abs(ah[~princ])) <= math.sqrt(q) / ph * (1 + 1e-9)
check("A4 chi_q(4x - c) = sum_psi ahat(psi) psi(x) on units (q odd squarefree, (q, c) = 1)", ok_rec)
check("A5 |ahat(psi_0)| = 1/phi(q) and |ahat(psi)| <= sqrt(q)/phi(q) (15 moduli, r = 1,2,4,5,7)", ok_size)
ok = True
for l in [3, 5, 7, 11, 13]:
    ah, princ, rec = coeffs(l, 0)
    chars, prim, pr, units = char_table(l)
    leg = np.array([gmpy2.legendre(a, l) for a in range(l)], dtype=float)
    k = [i for i in range(len(chars)) if np.allclose(chars[i][units], leg[units])][0]
    ok &= abs(ah[k] - 1) < 1e-12 and np.sum(np.abs(ah)) - 1 < 1e-9
check("A6 l | c: ahat(psi) = [psi = Legendre] (the case l | r and r = 0)", ok)

# ---------------------------------------------------------------- B: main term
print("== B: main term")


def E_closed(l, r, N):
    if l == N:
        return 1.0
    if l == 2:
        return 2 / 3 if r % 2 else 4 / 3
    if r % l == 0:
        return l * l / (l * l - 1)
    return 1 - 1 / ((l - 1) ** 2 * (l + 1))


def local_table(l, r, N, I):
    """e_l(i, j) for i <= I, j in {0, odd, even>=2}: means over p in (Z_l)^x of T_l(i, j; p)."""
    out = {}
    for i in range(I + 1):
        mod = l ** (2 * i + 3) if l == 2 else l ** (2 * i + 1)
        ps = np.array([p for p in range(mod) if p % l], dtype=np.int64)
        D = N * (4 * ps - r * r * N)
        ok = (D % (l ** (2 * i)) == 0)
        Dl = -(D // (l ** (2 * i)))
        if l == 2:
            ok &= (Dl % 4 == 0) | (Dl % 4 == 1)
            k1 = np.where(Dl % 2 == 0, 0, np.where((Dl % 8 == 1) | (Dl % 8 == 7), 1, -1))
        else:
            k1 = jacobi_vec(Dl, np.full(len(Dl), l))
        e0 = ok.mean(); eo = (ok * k1).mean(); ee = (ok * k1 * k1).mean()
        out[i] = (e0, eo, ee)
    return out


def local_sum(tab, l, principal_only=False):
    so = (1 / l) / (1 - 1 / l ** 2)      # sum_{j odd} l^-j
    se = (1 / l ** 2) / (1 - 1 / l ** 2)  # sum_{j even >= 2} l^-j
    s = 0.0
    for i, (e0, eo, ee) in tab.items():
        s += l ** (-i) * (e0 + (0 if principal_only else eo * so) + ee * se)
    return s


N = N0
ok = True; worst = 0.0
for r in range(0, 7):
    for l in [2, 3, 5, 7, 11, 13]:
        I = {2: 6, 3: 4, 5: 3, 7: 2}.get(l, 1)
        tab = local_table(l, r, N, I)
        F = local_sum(tab, l)
        tail = 32.0 * l ** (1 - 3 * (I + 1))  # sum_{i>I} l^-i * P(l^2i | D) * (1 + 1/(l-1)), P <= 8 l^(1-2i)
        worst = max(worst, abs(F - E_closed(l, r, N)))
        ok &= abs(F - E_closed(l, r, N)) <= tail + 1e-12
check("B1 local factor sum_{i,j} e_l(i,j) l^(-i-j) = E_l(r) (T5 closed forms), l <= 13, r = 0..6",
      ok, f"max |diff| {worst:.1e} (rigorous i-tail bound)")

# global double sum over m <= 40, n <= 20000 vs E(r) and nu'(r); principal-only variant
NMAX = 20000; MMAX = 40
spf = np.zeros(NMAX + 1, dtype=np.int64)
for p in sieve_primes(NMAX):
    blk = spf[p::p]; blk[blk == 0] = p
small_l = [int(l) for l in sieve_primes(50)]
big_primes = [int(l) for l in sieve_primes(10 ** 6)]


def fac_small(n):
    f = {}
    while n > 1:
        p = int(spf[n]); f[p] = f.get(p, 0) + 1; n //= p
    return f


res = {}
for r in range(0, 7):
    tabs = {l: local_table(l, r, N, {2: 6, 3: 4, 5: 3, 7: 2}.get(l, 1)) for l in small_l}

    def e_loc(l, i, j):
        if l in tabs:
            e0, eo, ee = tabs[l][i]
        else:  # l > 50, i = 0 (m <= 40); closed forms (checked against the tables for l <= 50 below)
            if i > 0:
                raise ValueError
            if l == N:
                e0, eo, ee = 1.0, 0.0, 0.0
            elif r % l == 0:
                e0, eo, ee = 1.0, 0.0, 1.0
            else:
                e0, eo, ee = 1.0, -1.0 / (l - 1), 1 - 1.0 / (l - 1)
        return e0 if j == 0 else (eo if j % 2 else ee)

    full = 0.0; princ = 0.0
    for m in range(1, MMAX + 1):
        fm = fac_small(m)
        em = 1.0
        for l, i in fm.items():
            em *= e_loc(l, i, 0)
        if em == 0:
            continue
        for n in range(1, NMAX + 1):
            fn = fac_small(n)
            e = em
            pr = True
            for l in set(fm) | set(fn):
                j = fn.get(l, 0); i = fm.get(l, 0)
                if j % 2:
                    pr = False
                if j:
                    e = e / e_loc(l, i, 0) * e_loc(l, i, j) if e_loc(l, i, 0) != 0 else 0.0
                if e == 0:
                    break
            if e:
                full += e / (m * n)
                if pr:
                    princ += e / (m * n)
    Er = 1.0
    for l in big_primes:
        Er *= E_closed(l, r, N)
    res[r] = (full, princ, Er)
# closed forms for l > 50 against enumeration at l = 53, 59, 61
ok = True
for l in [53, 59, 61]:
    for r in [0, 1, 2, 3]:
        e0, eo, ee = local_table(l, r, N, 0)[0]
        want = (1.0, 0.0, 1.0) if r % l == 0 else (1.0, -1 / (l - 1), 1 - 1 / (l - 1))
        ok &= max(abs(e0 - want[0]), abs(eo - want[1]), abs(ee - want[2])) < 1e-12
check("B2 closed forms e_l(0, j) = 1, -1/(l-1), 1 - 1/(l-1) (l not dividing 2rN), 1, 0, 1 (l | r)", ok)
tol = 3.0 / math.sqrt(NMAX)
for r in range(0, 7):
    full, princ, Er = res[r]
    nup = 1.0 if r == 0 else Er * (1 - 1 / ((N - 1) ** 2 * (N + 1)))
    ref = Er
    check(f"B3 r={r}: sum_(m<=40, n<=2e4) e(m,n)/(mn) = E(r)", abs(full - ref) < tol,
          f"sum {full:.5f}, E(r) {ref:.5f}" + (f" = zeta(2)(1-N^-2) {math.pi**2/6*(1-N**-2):.5f}" if r == 0 else
                                                f", nu'(r) {nup:.5f}") + f"; principal n only {princ:.5f}")
for r in range(1, 7):
    full, princ, Er = res[r]
    check(f"B4 r={r}: principal-only sum differs from E(r) (non-principal n carry main term)",
          abs(princ - Er) > 3 * tol, f"diff {princ - Er:+.4f}, tolerance {tol:.4f}")
full, princ, Er = res[0]
check("B5 r=0: principal-only sum = full sum (ahat(psi_0) = 0 for c = 0)", abs(princ - full) < 1e-12)

# ---------------------------------------------------------------- C: truncation per discriminant, N = 10007
print("== C: truncation of the Dirichlet series, N = 10007")
primes_all = sieve_primes(2 * N + 10)
JP = primes_all[(primes_all > N) & (primes_all <= 2 * N)].astype(np.int64)
logp = np.log(JP.astype(float))
theta = logp.sum()
x_len = N  # length of J
Ms = {"N^0.5": N ** 0.5, "N^0.75": N ** 0.75, "N^(7/8)": N ** 0.875, "N^1": N ** 1.0, "N^1.1": N ** 1.1}
MM = int(N ** 1.1) + 1
narr = np.arange(1, MM + 1, dtype=np.int64)
samp = JP[:: max(1, len(JP) // 40)][:40]
worst = {k: 0.0 for k in Ms}; ratioPV = 0.0; ratioB = 0.0; sqm = {k: 0.0 for k in Ms}
for r in [0, 1, 2]:
    for p in samp:
        p = int(p)
        D = N * (4 * p - r * r * N)
        Lstar = math.pi * float(pari.qfbhclassno(D)) / math.sqrt(D)
        parts = []
        for m in disc_ms(D):
            Dl = -D // (m * m)
            kv = kron_vec(Dl, narr).astype(float) / narr
            cs = np.cumsum(kv)
            D0 = int(pari.coredisc(Dl)); f = int(math.isqrt(Dl // D0))
            parts.append((m, cs, abs(D0), len(factor(f)) if f > 1 else 0))
        for k, M in Ms.items():
            Mi = int(M)
            S = sum(cs[Mi - 1] / m for (m, cs, q0, w) in parts)
            err = abs(S - Lstar)
            worst[k] = max(worst[k], err)
            pv = sum(2 * 2 ** w * math.sqrt(q0) * math.log(q0) / (m * Mi) for (m, cs, q0, w) in parts)
            bu = sum(2 ** w * Mi ** -0.5 * q0 ** (3 / 16) / m for (m, cs, q0, w) in parts)
            ratioPV = max(ratioPV, err / pv)
            ratioB = max(ratioB, err / bu)
            sqm[k] = max(sqm[k], err * math.sqrt(Mi))
print("   max |truncated - L*| over 120 (p, r):", {k: f"{v:.2e}" for k, v in worst.items()})
ks = ["N^0.75", "N^(7/8)", "N^1", "N^1.1"]
check("C1 max truncation error * sqrt(M) is non-increasing for M = N^0.75, N^(7/8), N, N^1.1 (square-root decay or better)",
      all(sqm[ks[i + 1]] <= sqm[ks[i]] for i in range(3)), " ".join(f"{k}:{sqm[k]:.2f}" for k in ks) +
      f"; max ratio to the Polya-Vinogradov shape {ratioPV:.4f} (that bound is > 1, i.e. vacuous, at this size)")
check("C2 truncation error <= Burgess r=2 shape 2^w(f) M^(-1/2) |D0|^(3/16) (constant 1, empirical)", ratioB <= 1,
      f"max ratio {ratioB:.3f}")

# mean of L*(D_r(p)) over all p in J against E(r); the same for the series truncated at M = N^(7/8) (step U2)
Mi = int(N ** 0.875); narr2 = np.arange(1, Mi + 1, dtype=np.int64)
for r in [0, 1, 2]:
    vals = np.array([math.pi * float(pari.qfbhclassno(N * (4 * int(p) - r * r * N))) / math.sqrt(N * (4 * int(p) - r * r * N))
                     for p in JP])
    tr = []
    for p in JP:
        D = N * (4 * int(p) - r * r * N)
        tr.append(sum(float((kron_vec(-D // (m * m), narr2) / narr2).sum()) / m for m in disc_ms(D)))
    tr = np.array(tr)
    mt = float((logp * tr).sum() / theta); sdt = float(tr.std()) / math.sqrt(len(tr))
    zt = (mt - res[r][2]) / sdt
    check(f"C3 r={r}: log-weighted mean over all {len(JP)} p in (N,2N] of the series truncated at M = N^(7/8) vs E(r)",
          abs(zt) < 3, f"mean {mt:.4f}, E(r) {res[r][2]:.4f}, z {zt:+.2f}; mean |truncation error| {np.abs(tr - vals).mean():.4f}")
    mean = float((logp * vals).sum() / theta)
    sd = float(vals.std()) / math.sqrt(len(vals))
    Er = res[r][2]
    z = (mean - Er) / sd
    check(f"C4 r={r}: log-weighted mean of L*(D_r(p)) over {len(JP)} primes in (N,2N] vs E(r)", abs(z) < 3,
          f"mean {mean:.4f}, E(r) {Er:.4f}, z {z:+.2f} (naive sd)")

# ---------------------------------------------------------------- D: ranges (ii)-(iv)
print("== D: shifted-prime character sums, N = 10007, J = (N, 2N]")


def sqfree_odd(lo, hi, avoid):
    out = []
    for q in range(lo, hi + 1):
        if q % 2 == 0 or q % avoid == 0:
            continue
        f = factor(q)
        if all(e == 1 for e in f.values()):
            out.append(q)
    return out


def S_shift(q, c):
    return float((logp * jacobi_vec(4 * JP - c, np.full(len(JP), q))).sum())


def a0(q, c):
    v = 1.0
    for l in factor(q):
        v *= -gmpy2.legendre(-c, l) / (l - 1) if c % l else 0.0
    return v


# (ii) small conductors: S_q = ahat(psi_0) theta + small
for r in [1, 2]:
    c = r * r * N
    devs = []
    for q in sqfree_odd(3, 31, N):
        if math.gcd(q, r) > 1:
            continue
        devs.append(abs(S_shift(q, c) - a0(q, c) * theta) / theta)
    check(f"D1 range (ii), r={r}: |S_q - ahat(psi0) theta|/theta small for odd squarefree q <= 31",
          max(devs) < 0.05, f"max {max(devs):.4f}; sqrt(log x / x) = {math.sqrt(math.log(N)/N):.4f}")

# large sieve inequality and the chain of step U5 on real data (Q = 16, 32, 64)
allchars = {}
for d in range(1, 129):
    allchars[d] = char_table(d)
Sp = {}
for d in range(1, 129):
    chars, prim, princ, units = allchars[d]
    vals = chars[:, JP % d]
    Sp[d] = (vals * logp).sum(axis=1)
lhs = sum(d / phi(d) * float(np.sum(np.abs(Sp[d][allchars[d][1]]) ** 2)) for d in range(1, 129))
rhs = (x_len - 1 + 128 ** 2) * float((logp ** 2).sum())
check("D2 multiplicative large sieve on J: sum_(d<=128) d/phi(d) sum* |S(psi)|^2 <= (|J| - 1 + Q^2) sum log^2 p",
      lhs <= rhs, f"lhs/rhs = {lhs/rhs:.3f}")
for r in [1, 0]:
    c = r * r * N
    for Q in [16, 32, 64]:
        qs = [q for q in sqfree_odd(Q + 1, 2 * Q, N) if math.gcd(q, max(r, 1)) == 1 or r == 0]
        L1 = sum(abs(S_shift(q, c)) / q for q in qs)
        if r == 1:
            L2 = 0.0; L3 = 0.0
            for q in qs:
                ah, princ, rec = coeffs(q, c)
                L2 += float(np.sum(np.abs(ah) * np.abs(Sp[q]))) / q
                L3 += (theta / phi(q) + math.sqrt(q) / phi(q) * float(np.sum(np.abs(Sp[q][~princ])))) / q
            # group by the primitive character: w_d = sum_{q ~ Q, d | q} q^(-1/2)/phi(q); C-S; large sieve
            w = {}
            for q in qs:
                for d in range(3, q + 1):
                    if q % d == 0:
                        w[d] = w.get(d, 0.0) + q ** -0.5 / phi(q)
            mainpart = sum(theta / (q * phi(q)) for q in qs)
            L4 = mainpart + sum(w[d] * float(np.sum(np.abs(Sp[d][allchars[d][1]]))) for d in w)
            cs1 = sum(w[d] ** 2 * int(allchars[d][1].sum()) * phi(d) / d for d in w)
            cs2 = sum(d / phi(d) * float(np.sum(np.abs(Sp[d][allchars[d][1]]) ** 2)) for d in w)
            L5 = mainpart + math.sqrt(cs1 * cs2)
            L6 = mainpart + math.sqrt(cs1 * (x_len - 1 + (2 * Q) ** 2) * float((logp ** 2).sum()))
            okc = L1 <= L2 * (1 + 1e-9) and L2 <= L3 * (1 + 1e-9) and abs(L3 - L4) <= 1e-6 * L3 and L4 <= L5 * (1 + 1e-9) and L5 <= L6
            check(f"D3 r=1, Q={Q}: chain L1 <= L2 <= L3 = L4 <= L5 <= L6 of step U5", okc,
                  f"L1/theta {L1/theta:.4f}, L3 {L3/theta:.4f}, L5 {L5/theta:.4f}, L6 {L6/theta:.4f}; "
                  f"Q^-1/2+(Q/x)^1/2 = {Q**-0.5+(Q/x_len)**0.5:.3f}")
        else:
            # r = 0: S_q = S(chi_q), chi_q primitive real; Cauchy-Schwarz + large sieve
            cs = math.sqrt(sum(1 / q ** 2 for q in qs)) * math.sqrt((x_len - 1 + (2 * Q) ** 2) * float((logp ** 2).sum()))
            check(f"D3 r=0, Q={Q}: sum (1/q)|S(chi_q)| <= (sum q^-2)^1/2 ((|J|-1+4Q^2) sum log^2 p)^1/2", L1 <= cs,
                  f"L1/theta {L1/theta:.4f}, bound/theta {cs/theta:.4f}")

# (iii) up to M = N^(7/8): block sums against the shape theta (Q^-1/2 + (Q/x)^1/2)
Mtr = int(N ** 0.875)
for r in [0, 1, 2]:
    c = r * r * N
    Q = 32; ok = True; rows = []
    while Q < Mtr:
        qs = [q for q in sqfree_odd(Q + 1, min(2 * Q, Mtr), N) if r == 0 or math.gcd(q, r) == 1]
        A = sum(abs(S_shift(q, c)) / q for q in qs)
        bound = theta * (Q ** -0.5 + (Q / x_len) ** 0.5)
        ok &= A <= bound
        rows.append(f"{Q}:{A/theta:.4f}")
        Q *= 2
    check(f"D4 range (iii), r={r}: sum_(q~Q) |S_q|/q <= theta (Q^-1/2 + (Q/x)^1/2), 32 <= Q < N^(7/8)", ok,
          "A(Q)/theta " + " ".join(rows))

# (iv) conductors near and above N (not needed after the Burgess truncation): sizes only
for r in [0, 1]:
    c = r * r * N
    qs = sqfree_odd(N - 300, N + 300, N)[:150] + sqfree_odd(int(N ** 1.1) - 300, int(N ** 1.1) + 300, N)[:150]
    vals = [abs(S_shift(q, c)) for q in qs]
    rel = max(vals) / math.sqrt(theta * math.log(theta))
    check(f"D5 range (iv) (empirical, not used): r={r}, 300 conductors near N and N^1.1: max|S_q| / sqrt(theta log theta)",
          rel < 6, f"{rel:.2f}; max|S_q|/theta = {max(vals)/theta:.4f}")

print(f"run time {time.time() - T0:.0f} s")
print(f"{NCHK} checks, {NFAIL} failed")
