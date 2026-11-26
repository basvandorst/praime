"""Task 2(a): own recomputation of the round-156 shifted-prime sums for the Hecke character of Q(i) (lambda(n) = (1/4) sum over
a^2 + b^2 = n of Re((a + bi)^4)/n^2, i.e. 2 cos 4theta at split p) on primes p in [L0, L1):
   S_h = sum_p lambda(p + h) ln p (p + h)^{-1/4} Hann(sqrt(p + h); V0, V1) e(-2 sqrt((p + h)/4)),  h = 2, 6,
lambda by lattice-point enumeration (no factoring), primes by an own segmented sieve.  Compared chunk by chunk (1e7) with
sums_holdout.jsonl.gz, and turned into a reading S / (C MT_sub) with the main term integrated over the sub-range only.
usage: t2_shift.py L0 L1"""
import sys, math, json, gzip
import numpy as np
L0, L1 = int(float(sys.argv[1])), int(float(sys.argv[2]))
V0, V1 = math.sqrt(3e7), math.sqrt(1e9); NN = 4; SEG = 5 * 10 ** 6
small = np.ones(26500, bool); small[:2] = False
for q in range(2, 163):
    if small[q]: small[q * q::q] = False
SP = np.nonzero(small)[0]

def primes_in(L, U):
    s = np.ones(U - L, bool)
    for q in SP:
        q = int(q)
        if q * q >= U: break
        st = max(q * q, ((L + q - 1) // q) * q); s[st - L::q] = False
    return np.nonzero(s)[0].astype(np.int64) + L

def isqrt_arr(x):
    r = np.floor(np.sqrt(x.astype(np.float64))).astype(np.int64)
    r[r * r > x] -= 1; r[(r + 1) * (r + 1) <= x] += 1
    return r

def lam_seg(L, U):
    """lambda(n) for n in [L, U), first-quadrant points with multiplicities m(a) m(b), m(0) = 1, m(x > 0) = 2."""
    amax = math.isqrt(U - 1); a = np.arange(amax + 1, dtype=np.int64)
    hi = isqrt_arr(U - 1 - a * a)
    lo2 = L - a * a; lo = np.where(lo2 <= 0, 0, isqrt_arr(np.maximum(lo2 - 1, 0)) + 1)
    cnt = np.maximum(hi - lo + 1, 0); tot = int(cnt.sum())
    A = np.repeat(a, cnt); off = np.repeat(np.cumsum(cnt) - cnt, cnt)
    B = np.repeat(lo, cnt) + (np.arange(tot, dtype=np.int64) - off); del off
    n = A * A + B * B
    af = A.astype(np.float64); bf = B.astype(np.float64)
    v = ((af * af - bf * bf) ** 2 - 4 * af * af * bf * bf) * np.where(A > 0, 2.0, 1.0) * np.where(B > 0, 2.0, 1.0)
    del A, B, af, bf
    out = np.bincount(n - L, weights=v, minlength=U - L)
    nn = np.arange(L, U, dtype=np.float64)
    return out / (4 * nn * nn)

def hann(v):
    u = (v - V0) / (V1 - V0); return np.where((u > 0) & (u < 1), np.sin(np.pi * u) ** 2, 0.0)

res = {}   # chunk (1e7) -> dict h -> [S, noise, count]
for L in range(L0, L1, SEG):
    U = L + SEG; pr = primes_in(L, U); lam = lam_seg(L, U + 8)
    key = (L // 10 ** 7) * 10 ** 7
    for h in (2, 6):
        n = pr + h; nf = n.astype(np.float64); w = lam[n - L] * np.log(pr.astype(np.float64)) * nf ** -0.25 * hann(np.sqrt(nf))
        s = complex(np.dot(w, np.exp(-2j * math.pi * 2 * np.sqrt(nf / NN))))
        r = res.setdefault(key, {}).setdefault(h, [0j, 0.0, 0]); r[0] += s; r[1] += float(w @ w); r[2] += len(n)
    if L == L0:   # sanity: lambda at small split primes in this segment equals 2 cos 4 theta (checked separately in t2_euler.py)
        print('segment', L, 'primes', len(pr), 'mean lambda(p+2)^2 %.4f' % float(np.mean(lam[pr + 2 - L] ** 2)), flush=True)

stored = {}
with gzip.open('/home/bas/praime/writeup/data/study/r156/sums_holdout.jsonl.gz', 'rt') as f:
    for line in f:
        d = json.loads(line)
        if L0 <= d['L'] < L1: stored[d['L']] = d['res']['A4']
tot = {2: [0j, 0.0, 0], 6: [0j, 0.0, 0]}; maxrel = {2: 0.0, 6: 0.0}
for key in sorted(res):
    for h in (2, 6):
        S, nz, c = res[key][h]; st = stored[key]['shift%d' % h]; Ss = complex(*st['smooth'][0])
        maxrel[h] = max(maxrel[h], abs(S - Ss) / math.sqrt(nz / 2)); assert c == st['count'], (key, h, c, st['count'])
        tot[h][0] += S; tot[h][1] += nz; tot[h][2] += c
print('chunks compared:', len(res), '; counts identical; largest |own - stored| in noise units: h=2 %.1e, h=6 %.1e' % (maxrel[2], maxrel[6]))
# main term over the sub-range: (3/2) K1 int W(v) (1 - i kappa sqrt(N) / (2 pi v)) dv, K1 = (2/3) e^{-i pi/4}, kappa = -63/16
K1 = (2 / 3) * complex(math.cos(math.pi / 4), -math.sin(math.pi / 4)); kappa = 3 / 48 + 0.5 * ((2 - 4) + (3 - 9))
for h in (2, 6):
    v = np.linspace(math.sqrt(L0 + h), math.sqrt(L1 + h), 200001)
    f = hann(v) * (1 - 1j * kappa * math.sqrt(NN) / (2 * math.pi * v)); MT = 1.5 * K1 * np.trapezoid(f, v)
    print('h=%d  sub-range [%.0e, %.0e): %d primes  S = %.2f %+.2fi  MT_sub = %.2f %+.2fi  S/MT = %.4f %+.4fi  noise error of S/MT %.4f'
          % (h, L0, L1, tot[h][2], tot[h][0].real, tot[h][0].imag, MT.real, MT.imag, (tot[h][0] / MT).real, (tot[h][0] / MT).imag,
             math.sqrt(tot[h][1] / 2) / abs(MT)))
json.dump({str(h): [tot[h][0].real, tot[h][0].imag, tot[h][1], tot[h][2]] for h in (2, 6)}, open('t2_shift_out.json', 'w'))
