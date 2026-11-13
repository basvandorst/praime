"""The hidden part of the class sums: every residue class mod q carries its whole progression.

C(n) = (1/(q-1)) sum_chi chibar(n) log L(2, chi) = sum_{j >= 0} a(n + jq) / (n + jq)^2,  a(m) = Lambda(m)/ln m
(a = 1 for a prime, 1/k for p^k, 0 otherwise). Peeling: for j = 0, 1, 2, ... take r = C(n) - (terms already
decided), estimate a(n + jq) ~ r (n + jq)^2, decide 'prime' if the estimate exceeds 3/4, subtract the decided term
(prime powers p^k get 1/k from the decision on p, which is always made earlier since p <= sqrt(n + jq) < q).
Not a proof beyond n < q/3: the rest of the class can in principle add up to more than 1/4.
Usage: charlog_peel.py q J
"""
import sys, math, time
import numpy as np
import scipy.fft as sfft
from scipy.special import zeta
from sympy import isprime, primitive_root
import gmpy2

q = int(sys.argv[1]); J = int(sys.argv[2]); s = 2.0
assert isprime(q)
t0 = time.time()
g = int(primitive_root(q))
B = 1 << 16
pw = np.empty(q - 1, dtype=np.int32); v = 1
for j in range(min(B, q - 1)):
    pw[j] = v; v = v * g % q
gB = pow(g, B, q)
for st in range(B, q - 1, B):
    en = min(st + B, q - 1); pw[st:en] = (pw[st - B:en - B].astype(np.int64) * gB) % q
H = zeta(s, pw / q)
L = sfft.ifft(H) * (q - 1) * float(q) ** (-s); del H
np.log(L, out=L)
Cg = sfft.fft(L, overwrite_x=True) / (q - 1); del L
C = np.zeros(q); C[pw] = Cg.real                    # C[n] for every class n = 1 .. q-1
imag = np.zeros(q); imag[pw] = Cg.imag; del Cg, pw
t1 = time.time()
M = (J + 1) * q
truth = np.fromiter((gmpy2.is_prime(v) for v in range(M)), dtype=bool, count=M)   # reference only, each number on its own
t2 = time.time()
# peel all classes at once, band by band
dec = np.zeros(M, dtype=bool)
rem = C.copy()
PP = sorted({(r ** k, r, k) for r in range(2, int(math.isqrt(M)) + 1) for k in range(2, 64) if r ** k < M and not any(gmpy2.iroot(r, e)[1] for e in range(2, 32))})
lines = [f'q = {q}: L-values, logs, FFTs {t1 - t0:.1f} s; reference primes below {M:,} {t2 - t1:.0f} s']
CH = 5_000_000
wrong_all = {}
for j in range(J + 1):
    tp = dp = fp = mp_ = 0; first = []; wl = []
    ppj = [(val, root, kk) for (val, root, kk) in PP if j * q < val < (j + 1) * q]
    for lo in range(1, q, CH):
        nn = np.arange(lo, min(lo + CH, q), dtype=np.int64)
        m = nn + j * q; mf = m.astype(float) ** 2
        est = rem[nn] * mf
        if j == 0 and lo == 1:
            lines.append(f'largest |Im C| (n+q)^2 {float(np.max(np.abs(imag[1:]) * (np.arange(1, q) + q) ** 2.0)):.2e}')
        pr = (est > 0.75) & ((m % 2 == 1) | (m == 2))
        loc = [(val - j * q - lo, root, kk) for (val, root, kk) in ppj if lo <= val - j * q < lo + len(nn)]
        for (i_, root, kk) in loc: pr[i_] = False
        dec[m] = pr
        a_d = pr.astype(float)
        for (i_, root, kk) in loc: a_d[i_] = (1.0 / kk) if dec[root] else 0.0
        rem[nn] -= a_d / mf
        t = truth[m]; w = np.nonzero(pr != t)[0]
        tp += int(t.sum()); dp += int(pr.sum()); fp += int(np.sum(pr & ~t)); mp_ += int(np.sum(~pr & t))
        if len(first) < 5: first += [int(x) for x in m[w[:5 - len(first)]]]
        if len(wl) < 20000: wl += [int(x) for x in nn[w[:20000 - len(wl)]]]
    wrong_all[j] = np.array(wl, dtype=np.int64)
    lines.append(f'band {j}: numbers {j * q + 1:,} .. {(j + 1) * q - 1:,}: {tp:,} primes; decided prime {dp:,}; wrong {fp + mp_} '
                 f'(false prime {fp}, missed prime {mp_}); first wrong {first}')
    print(lines[-1], flush=True)
if 'imag' in dir(): del imag
# every miss: is it explained by later primes in the same class?  rest = sum_{j' > j} a(n + j'q) ((n + jq)/(n + j'q))^2
expl = []
n = None
for j in range(1, J + 1):
    nw = wrong_all[j]
    if len(nw) == 0: continue
    wrong = nw
    rest = np.zeros(len(wrong)); later = []
    for jj in range(j + 1, J + 40):
        mm = (wrong + jj * q)
        aa = np.array([1.0 if gmpy2.is_prime(int(x)) else 0.0 for x in mm])
        rest += aa * ((wrong + j * q) / mm.astype(float)) ** 2
        if jj <= j + 4: later.append(aa)
    later = np.array(later).T
    pat = {}
    for row in later.astype(int):
        key = ''.join(map(str, row)); pat[key] = pat.get(key, 0) + 1
    top = sorted(pat.items(), key=lambda kv: -kv[1])[:4]
    expl.append(f'band {j}: {len(wrong)} misses; rest of the class from the true primes: min {rest.min():.3f}, median {np.median(rest):.3f} '
                f'(misses with rest > 0.75 - 0.01: {int(np.sum(rest > 0.74))}); which of n+(j+1)q .. n+(j+4)q are prime, most common: {top}')
lines += expl
np.save(f'peel_{q}_dec.npy', np.packbits(dec))
print('\n'.join(lines))
open(f'peel_{q}.txt', 'w').write('\n'.join(lines) + '\n')
