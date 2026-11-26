"""Task 1(b): trivial comparators for round 157 on the same 40 curves.
(A) local constraints only: bad primes must have a(l) in {0, +-1}; exponent 1 if |a(l)| = 1, exactly 2 if a(l) = 0 and l >= 5
    (2..8 at l = 2, 2..5 at l = 3).  Count the candidates N' <= 6000 that survive.
(B) functional-equation check (the classical method, Dokchitser 2004 / lfuncheckfeq): for every N' in [2, 6000], eps = +-1,
    build a_n, n <= M, from a(p), p <= M (Hecke at p !| N', a(p^k) = a(p)^k at p | N'), phi(t) = sum a_n exp(-2 pi n t / sqrt N'),
    discrepancy D = max_t |phi(1/t) - eps t^2 phi(t)| / max(|phi(1/t)|, t^2 |phi(t)|), t in {1.15, 1.3, 1.45}.  Estimate = argmin D.
a(p) by own point counting (numpy); the truth (N, eps) from PARI only for scoring; PARI ellap compared with the own counts.
"""
import json, math, sys
import numpy as np
import cypari2
pari = cypari2.Pari()
CUR = json.load(open('/home/bas/praime/writeup/code/study/resonance/spectro/testset157.json'))
MMAX = 1000; NMAX = 6000

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for q in range(2, int(n ** 0.5) + 1):
        if s[q]: s[q * q::q] = False
    return np.nonzero(s)[0]
PR = primes_upto(NMAX)

def ap_count(ai, p):
    a1, a2, a3, a4, a6 = ai
    x = np.arange(p, dtype=np.int64)
    if p == 2:
        cnt = 0
        for xx in range(2):
            for yy in range(2):
                if (yy * yy + a1 * xx * yy + a3 * yy - (xx ** 3 + a2 * xx * xx + a4 * xx + a6)) % 2 == 0: cnt += 1
        return 2 - cnt
    D = ((a1 * x + a3) ** 2 + 4 * (x ** 3 + a2 * x * x + a4 * x + a6)) % p
    leg = -np.ones(p, np.int64); leg[(np.arange(1, p) ** 2) % p] = 1; leg[0] = 0
    return int(-leg[D].sum())

def spf_table(n):
    s = np.zeros(n + 1, np.int64)
    for q in range(2, n + 1):
        if s[q] == 0: s[q::q][s[q::q] == 0] = q
    return s

def run(M, report=False):
    spf = spf_table(M)
    n = np.arange(M + 1)
    # index sets per prime l <= M: (k, n = l^k m, m) with gcd(m, l) = 1
    idx = {}
    for l in PR[PR <= M]:
        l = int(l); lst = []; k = 1; lk = l
        while lk <= M:
            m = np.arange(1, M // lk + 1); m = m[m % l != 0]; lst.append((k, lk * m, m)); k += 1; lk *= l
        idx[l] = lst
    ts = np.array([1.15, 1.3, 1.45]); res = {}
    for name, ai in CUR.items():
        ap = {int(p): ap_count(ai, int(p)) for p in PR[PR <= max(M, 97)]}
        # all-good a_n
        A0 = np.zeros(M + 1); A0[1] = 1.0
        for m in range(2, M + 1):
            q = int(spf[m]); k = 0; r = m
            while r % q == 0: r //= q; k += 1
            # a(q^k) good
            t0, t1 = 1.0, float(ap[q])
            for _ in range(k - 1): t0, t1 = t1, ap[q] * t1 - q * t0
            A0[m] = A0[r] * (t1 if k >= 1 else 1.0)
        best = []
        for Np in range(2, NMAX + 1):
            bad = [l for l in PR[:200] if Np % l == 0] if Np > 1 else []
            r = Np
            for l in bad:
                while r % l == 0: r //= l
            if r > 1: bad.append(r)
            A = A0
            badM = [l for l in bad if l <= M]
            if badM:
                A = A0.copy()
                for l in badM:
                    for k, nk, m in idx[l]: A[nk] = A[m] * float(ap[l]) ** k
            sq = math.sqrt(Np)
            Et = np.exp(-2 * math.pi * np.outer(ts, n) / sq); Ei = np.exp(-2 * math.pi * np.outer(1 / ts, n) / sq)
            ph_t = Et @ A; ph_i = Ei @ A
            for eps in (1, -1):
                D = np.max(np.abs(ph_i - eps * ts ** 2 * ph_t) / np.maximum(np.abs(ph_i), ts ** 2 * np.abs(ph_t)))
                best.append((D, Np, eps))
        best.sort()
        E = pari.ellinit(ai); Nt = int(pari.ellglobalred(E)[0]); et = int(pari.ellrootno(E))
        rank = [i for i, b in enumerate(best) if b[1] == Nt and b[2] == et][0]
        Dt = best[rank][0]; Dnext = [b[0] for b in best if not (b[1] == Nt and b[2] == et)][0]
        res[name] = (Nt, et, best[0][1], best[0][2], rank, Dt, Dnext)
        if report:
            # (A) local candidate count and PARI comparison
            apall = {int(p): ap_count(ai, int(p)) for p in PR}
            B = {l: v for l, v in apall.items() if abs(v) <= 1}
            ok = []
            for Np in range(2, NMAX + 1):
                r = Np; good = True
                for l in PR:
                    l = int(l)
                    if l > r: break
                    if r % l: continue
                    e = 0
                    while r % l == 0: r //= l; e += 1
                    if l not in B: good = False; break
                    if abs(B[l]) == 1 and e != 1: good = False; break
                    if B[l] == 0 and ((l >= 5 and e != 2) or (l == 2 and not 2 <= e <= 8) or (l == 3 and not 2 <= e <= 5)): good = False; break
                if good: ok.append(Np)
            pa = [int(pari.ellap(E, p)) for p in PR[PR <= M]]
            mism = sum(1 for p, x in zip(PR[PR <= M], pa) if ap[int(p)] != x)
            u = pari.ellglobalred(E)[1][0]
            res[name] += (len(ok), Nt in ok, ok.index(Nt) if Nt in ok else -1, mism, int(u))
    return res

if __name__ == '__main__':
    for M in (int(x) for x in sys.argv[1:]):
        rep = (M == MMAX)
        res = run(M, report=rep)
        hits = sum(1 for v in res.values() if v[0] == v[2] and v[1] == v[3])
        print('M = %d (a(p) for p <= %d): FE hits %d of 40; truth D max %.1e; next-best D min %.1e' %
              (M, M, hits, max(v[5] for v in res.values()), min(v[6] for v in res.values())))
        if rep:
            for k, v in sorted(res.items(), key=lambda kv: kv[1][0]):
                print('%-9s N=%4d eps=%+d FE->%4d%+d D_true %.0e D_next %.0e | local cands %4d truth in: %s (pos %d) | ellap mism %d u=%d'
                      % (k, v[0], v[1], v[2], v[3], v[5], v[6], v[7], v[8], v[9], v[10], v[11]))
            print('local candidates: median %d, min %d, max %d' % (np.median([v[7] for v in res.values()]), min(v[7] for v in res.values()), max(v[7] for v in res.values())))
