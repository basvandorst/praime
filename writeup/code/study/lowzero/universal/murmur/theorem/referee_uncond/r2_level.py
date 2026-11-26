"""Referee: the unconditional chain at one prime level N = 20011, all primes p in (N, 2N], r = 0, 1, 2.
exact class-number mean  vs  truncated series (n <= N^(7/8))  vs  E(r);
decomposition of the truncated n-sum by the conductor q(n), by n'' and by the square/2-part cut (set script N_L);
empirical per-(m,n) means against the local prediction e(m,n). Own code."""
import math, time, sys
import numpy as np
import cypari2
from sympy import factorint as _fi, primerange
from functools import lru_cache
@lru_cache(maxsize=None)
def _fic(n):
    return dict(_fi(n))
def factorint(n):
    return _fic(n)

pari = cypari2.Pari()
T0 = time.time()
N = 20011
M = int(N ** 0.875)
primes = np.array(list(primerange(N + 1, 2 * N + 1)), dtype=np.int64)
logp = np.log(primes.astype(float)); theta = logp.sum()
small = [int(l) for l in primerange(2, M + 1)]
spf = np.zeros(M + 1, dtype=np.int64)
for l in small:
    blk = spf[l::l]; blk[blk == 0] = l
print(f"N={N}, M=N^(7/8)={M}, {len(primes)} primes in (N,2N], theta={theta:.1f}")

QR = {}
for l in small:
    if l > 2:
        t = np.full(l, -1, dtype=np.int8); t[0] = 0
        t[(np.arange(1, l) ** 2) % l] = 1
        QR[l] = t


def kron_matrix(Deltas):
    """K[n, j] = (Delta_j / n) for 1 <= n <= M (Kronecker), Deltas int64 array (all = 0,1 mod 4)."""
    Dl = np.asarray(Deltas, dtype=np.int64)
    K = np.zeros((M + 1, len(Dl)), dtype=np.int8)
    K[1] = 1
    for l in small:
        if l == 2:
            r8 = Dl % 8
            K[2] = np.where(Dl % 2 == 0, 0, np.where((r8 == 1) | (r8 == 7), 1, -1))
        else:
            K[l] = QR[l][Dl % l]
    for n in range(2, M + 1):
        l = spf[n]
        if l != n:
            K[n] = K[l] * K[n // l]
    return K


def fac(n):
    f = {}
    while n > 1:
        l = int(spf[n]); f[l] = f.get(l, 0) + 1; n //= l
    return f


def kr(D, n):
    return int(pari.kronecker(D, n))


cache = {}


def e_loc(l, i, jc, r):
    """mean over p in (Z/l^k)^x of 1[l^2i | D, disc]*kappa_l(-D/l^2i)^j, jc in {0,'o','e'}"""
    key = (l, i, jc, r)
    if key in cache:
        return cache[key]
    mod = l ** (2 * i + 1) if l != 2 else 2 ** (2 * i + 4)
    ps = np.array([p for p in range(mod) if p % l], dtype=np.int64)
    ps = ps + mod * ((r * r * N) // (4 * mod) + 2)
    D = N * (4 * ps - r * r * N)
    ok = D % (l ** (2 * i)) == 0
    Dl = -(D // (l ** (2 * i)))
    if l == 2:
        ok &= (Dl % 4 == 0) | (Dl % 4 == 1)
        kap = np.where(Dl % 2 == 0, 0, np.where((Dl % 8 == 1) | (Dl % 8 == 7), 1, -1))
    else:
        kap = QR[l][Dl % l] if l in QR else np.array([kr(int(x), l) for x in Dl])
    w = ok * (1 if jc == 0 else (kap if jc == 'o' else kap * kap))
    v = float(w.mean())
    cache[key] = v
    return v


def e_mn(m, n, r):
    fm = factorint(m) if m > 1 else {}
    fn = fac(n)
    v = 1.0
    for l in set(fm) | set(fn) | {2}:
        i = fm.get(l, 0); j = fn.get(l, 0)
        jc = 0 if j == 0 else ('o' if j % 2 else 'e')
        if l == N:
            v *= 0.0 if (i or j) else 1.0
        elif i == 0 and l != 2 and (r % l != 0 or r == 0 and False):
            # l not dividing 2 m r N: closed forms (r >= 1) ; r = 0 handled below
            if r == 0:
                v *= 1.0 if jc != 'o' else 0.0   # kappa = (-4pN/l), p uniform: mean 0 for odd j, 1 for even
            else:
                v *= 1.0 if jc == 0 else (-1.0 / (l - 1) if jc == 'o' else 1 - 1.0 / (l - 1))
        elif i == 0 and l != 2:  # l | r, l odd: kappa = (-4pN/l) = +-1 equally
            v *= 1.0 if jc != 'o' else 0.0
        else:
            v *= e_loc(l, i, jc, r)
        if v == 0:
            break
    return v


def split(m, n, r):
    """(q, s, n'') with q = prod{l : v_l(n) odd, l not | 2 m r N}, q-natural = q s^2, n'' = n / (q s^2)."""
    f = fac(n)
    badset = set(factorint(2 * m * N)) | (set(factorint(r)) if r >= 1 else set())
    q = 1; qn = 1
    for l, e in f.items():
        if e % 2 and l not in badset:
            q *= l; qn *= l ** e
    s = math.isqrt(qn // q)
    return q, s, n // qn


narr = np.arange(1, M + 1, dtype=float)
for r in [0, 1, 2]:
    t1 = time.time()
    u = 4 * primes - r * r * N
    D = N * u
    Lstar = np.array([math.pi * float(pari.qfbhclassno(int(d))) / math.sqrt(int(d)) for d in D])
    rows_p = []; rows_m = []; rows_D = []
    for k, d in enumerate(D):
        d = int(d)
        core = int(pari.core(d))
        s2 = math.isqrt(d // core)
        for m in [int(x) for x in pari.divisors(s2)]:
            if (-d // (m * m)) % 4 in (0, 1):
                rows_p.append(k); rows_m.append(m); rows_D.append(-d // (m * m))
    rows_p = np.array(rows_p); rows_m = np.array(rows_m); rows_D = np.array(rows_D, dtype=np.int64)
    K = kron_matrix(rows_D)
    trunc_rows = (K[1:].astype(float) / narr[:, None]).sum(axis=0)
    S_M = np.zeros(len(primes)); np.add.at(S_M, rows_p, trunc_rows / rows_m)
    err = S_M - Lstar
    ex_mean = float((logp * Lstar).sum() / theta)
    tr_mean = float((logp * S_M).sum() / theta)
    sd = float(Lstar.std() / math.sqrt(len(primes)))
    Er = 1.0
    for l in primerange(2, 2 * 10 ** 6):
        if l == N:
            continue
        Er *= (2 / 3 if r % 2 else 4 / 3) if l == 2 else (l * l / (l * l - 1) if r % l == 0 else 1 - 1 / ((l - 1) ** 2 * (l + 1)))
    print(f"\n== r={r}: rows (p,m) {len(rows_p)}, m values {sorted(set(rows_m.tolist()))[:12]}..., max m {rows_m.max()}")
    print(f"  exact mean <L*> = {ex_mean:.5f}   truncated mean = {tr_mean:.5f}   E(r) = {Er:.5f}   naive sd {sd:.5f}")
    print(f"  per-p truncation error: max {np.abs(err).max():.4f}, rms {math.sqrt((err**2).mean()):.4f}, mean {float((logp*err).sum()/theta):+.5f}")
    # Burgess shape for the worst p: M^-1/2 |D0|^(3/16)
    print(f"  Burgess shape M^-1/2 (4N^2)^(3/16) = {M**-0.5*(4*N*N)**(3/16):.3f};  Polya-Vinogradov shape sqrt(4N^2)log/M = {2*N*math.log(4*N*N)/M:.1f}")
    # per-(m,n) empirical means
    ms = sorted(set(rows_m.tolist()))
    tot_emp = 0.0; tot_pred = 0.0
    blocks = {}  # key -> [emp, pred, abs_emp]
    nsec = {}
    Lcut = {4: [0.0, 0.0], 16: [0.0, 0.0], 64: [0.0, 0.0]}
    for m in ms:
        sel = rows_m == m
        wts = logp[rows_p[sel]]
        mu = (K[1:, sel].astype(float) * wts[None, :]).sum(axis=1) / theta  # mu[n-1]
        for n in range(1, M + 1):
            emp = mu[n - 1] / (m * n)
            pred = e_mn(m, n, r) / (m * n) if m <= 40 else 0.0
            tot_emp += emp; tot_pred += pred
            q, s, n2 = split(m, n, r)
            if q == 1:
                key = "q=1"
            elif q <= 32:
                key = "q in (1,32]"
            else:
                key = f"q~{2 ** int(math.log2(q))}"
            b = blocks.setdefault(key, [0.0, 0.0, 0.0]); b[0] += emp; b[1] += pred; b[2] += abs(emp)
            kk = "n''<=8" if n2 <= 8 else ("n''<=64" if n2 <= 64 else ("n''<=512" if n2 <= 512 else "n''>512"))
            c = nsec.setdefault(kk, [0.0, 0.0, 0.0]); c[0] += emp; c[1] += pred; c[2] += abs(emp)
            fo = fac(n); a2 = fo.get(2, 0)
            n2sq = 1
            for l, e in fo.items():
                if l != 2:
                    n2sq *= l ** (e // 2)
            for Lc in Lcut:
                if 2 ** a2 > Lc or n2sq > Lc:
                    Lcut[Lc][0] += emp; Lcut[Lc][1] += pred
    print(f"  sum over (m,n<=M) of empirical means /(mn) = {tot_emp:.5f} (= truncated mean);  predicted sum e(m,n)/(mn) (m<=40) = {tot_pred:.5f}")
    print("  by conductor block      empirical   predicted   sum|emp|")
    for key in sorted(blocks, key=lambda k: (0 if k == 'q=1' else (1 if k.startswith('q in') else 2), int(k.split('~')[1]) if '~' in k else 0)):
        b = blocks[key]
        print(f"    {key:16s} {b[0]:+10.5f} {b[1]:+10.5f} {b[2]:10.5f}")
    print("  by n'' size            empirical   predicted   sum|emp|")
    for key in ["n''<=8", "n''<=64", "n''<=512", "n''>512"]:
        if key in nsec:
            c = nsec[key]
            print(f"    {key:16s} {c[0]:+10.5f} {c[1]:+10.5f} {c[2]:10.5f}")
    for Lc, (a, b) in Lcut.items():
        print(f"  n outside N_L (2-part > L or odd square part > L), L={Lc}: empirical {a:+.5f}, predicted {b:+.5f}, "
              f"trivial per-p bound log(M)/L-ish {2*(math.log(M)+1)/Lc:.3f}")
    print(f"  time {time.time()-t1:.0f} s")
print(f"total {time.time()-T0:.0f} s")
