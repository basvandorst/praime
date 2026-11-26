# Referee: the unconditional chain at an even level, own code.
# For p prime in (N, 2N], r = 0, 1, 2: exact L*(D) = pi H(D)/sqrt(D) (PARI Hurwitz class number),
# the truncated series sum_{m in M(D)} m^-1 sum_{n <= N^(7/8)} xi_{-D/m^2}(n)/n (own Kronecker),
# split by the conductor q(n) = prod{l odd : v_l(n) odd, l not | m r N}; local predictions per block.
import sys, math
import numpy as np
from sympy import primerange, factorint, isprime
import cypari2
pari = cypari2.Pari()
pari.allocatemem(10**8)

Nlist = [int(a) for a in sys.argv[1:]] or [2 * 10007]

def kron2(a):
    if a % 2 == 0: return 0
    return 1 if a % 8 in (1, 7) else -1

for N in Nlist:
    PN = set(factorint(N))
    assert all(e == 1 for e in factorint(N).values())
    M = int(N ** 0.875)
    primesM = list(primerange(2, M + 1))
    # q0(n): product of odd primes with odd exponent
    q0 = np.ones(M + 1, dtype=np.int64)
    for l in primesM:
        if l == 2: continue
        k = 1; lk = l
        while lk <= M:
            # toggle: multiply by l when exponent becomes odd, divide when even
            idx = np.arange(lk, M + 1, lk)
            if k % 2: q0[idx] *= l
            else: q0[idx] //= l
            k += 1; lk *= l
    inv_n = np.zeros(M + 1); inv_n[1:] = 1.0 / np.arange(1, M + 1)
    qcache = {}
    def qarr(m, r):
        key = (m, r)
        if key not in qcache:
            excl = 1
            for l in set(factorint(m * max(r, 1) * N)) | PN:
                if l % 2: excl *= l
            g = np.gcd(q0, excl)
            qcache[key] = q0 // g
        return qcache[key]
    blocks = [(1, 1), (2, 32), (33, 1024), (1025, 10 ** 9)]

    def chi_arr(Delta):
        c = np.ones(M + 1, dtype=np.int8); c[0] = 0
        for l in primesM:
            v = kron2(Delta) if l == 2 else (0 if Delta % l == 0 else (1 if pow(Delta % l, (l - 1) // 2, l) == 1 else -1))
            if v == 1: continue
            if v == 0:
                c[l::l] = 0; continue
            lk = l
            while lk <= M:
                c[lk::lk] *= -1
                lk *= l
        return c

    ps = [p for p in primerange(N + 1, 2 * N + 1)]
    wts = np.array([math.log(p) for p in ps])
    print("N = %d = %s, M = N^(7/8) = %d, %d primes in (N, 2N]" % (N, sorted(PN), M, len(ps)))

    # local predictions
    def E_loc(l, r):
        if N % l == 0: return 1.0
        if r % l == 0: return l * l / (l * l - 1.0)
        return 1 - 1.0 / ((l - 1) ** 2 * (l + 1))
    LP = list(primerange(2, 2 * 10 ** 6))
    for r in [0, 1, 2]:
        E = 1.0; a_prod = 1.0
        ratio = {}
        for l in LP:
            El = E_loc(l, r); E *= El
            if l % 2 and N % l and (r == 0 or r % l):
                # q-part: i = 0, j odd: g_l = sum_{j odd} e_l(0,j) l^-j, e_l(0, odd) = -1/(l-1) for l not | 2rN
                # for r = 0: e_l(0, j odd) = mean of (-D/l)^j = mean of (-pN/l) = 0
                gl = (-l / ((l - 1.0) * (l * l - 1.0))) if r else 0.0
                al = El - gl
                a_prod *= al
                if l <= 1024: ratio[l] = gl / al
            else:
                a_prod *= El
        # sum over squarefree odd q <= 1024 coprime to rN of prod ratio
        V = {1: 1.0}
        for l in sorted(ratio):
            new = {}
            for q, v in V.items():
                if q * l <= 1024: new[q * l] = v * ratio[l]
            V.update(new)
        pred = [a_prod * sum(v for q, v in V.items() if lo <= q <= hi) for lo, hi in blocks[:3]]
        pred.append(E - sum(pred))
        # data
        Lstar = np.zeros(len(ps)); ser = np.zeros(len(ps)); blk = np.zeros((len(ps), 4))
        for t, p in enumerate(ps):
            D = N * (4 * p - r * r * N)
            H = pari.qfbhclassno(D)
            Lstar[t] = math.pi * float(H) / math.sqrt(D)
            f = factorint(D); base = 1
            for l, e in f.items(): base *= l ** (e // 2)
            ms = [m for m in range(1, base + 1) if base % m == 0 and D % (m * m) == 0 and (-D // (m * m)) % 4 in (0, 1)]
            assert all(math.gcd(m, N) == 1 for m in ms)
            for m in ms:
                c = chi_arr(-D // (m * m))
                if N % 2 == 0: assert not c[2::2].any()
                for q in PN: assert not c[q::q].any()
                w = c * inv_n / m
                qa = qarr(m, r)
                ser[t] += w.sum()
                for b, (lo, hi) in enumerate(blocks):
                    blk[t, b] += w[(qa >= lo) & (qa <= hi)].sum()
        def mean(x):
            mu = np.sum(wts * x) / wts.sum()
            se = np.std(x) / math.sqrt(len(x))
            return mu, se
        mL, sL = mean(Lstar); mS, sS = mean(ser)
        rms = math.sqrt(np.mean((ser - Lstar) ** 2))
        print(" r=%d  E^(N)(r) = %.5f | exact mean L* = %.5f (z %+.2f) | series cut at M: %.5f (z %+.2f), rms per-p cut error %.4f"
              % (r, E, mL, (mL - E) / sL, mS, (mS - E) / sS, rms))
        for b, (lo, hi) in enumerate(blocks):
            mb, sb = mean(blk[:, b])
            lab = "q = 1" if b == 0 else ("%d <= q <= %d" % (lo, hi) if hi < 10 ** 9 else "q > %d" % (lo - 1))
            print("     block %-16s data %+.5f (se %.5f)  predicted %+.5f%s" % (lab, mb, sb, pred[b], "  (= E - lower blocks, infinite tail)" if b == 3 else ""))
        sys.stdout.flush()
