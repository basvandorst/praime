"""Probe B: the least prime in every residue class modulo q = 30,000,001 (prime), and the greatest of them, P(q).
Each class a = 1 .. q-1 is walked along a, a + q, a + 2q, ...; candidates with a factor below 17 are skipped, the rest
are tested one by one with gmpy2 (BPSW; a proof below 2^64). Wagstaff (Math. Comp. 33, 1979) predicts P(q) about
phi(q) ln q ln phi(q)."""
import math, time, sys
import numpy as np
import gmpy2
q = int(sys.argv[1]); t0 = time.time()
assert gmpy2.is_prime(q)
rem = np.arange(1, q, dtype=np.int64)
jfound = np.zeros(q, dtype=np.int32)
small = (2, 3, 5, 7, 11, 13)
j = 0; hist = []
while len(rem):
    m = rem + j * q
    ok = np.ones(len(m), dtype=bool)
    for p in small:
        ok &= (m % p != 0) | (m == p)
    cand = np.nonzero(ok)[0]
    isp = np.fromiter((gmpy2.is_prime(int(v)) for v in m[cand]), dtype=bool, count=len(cand))
    hit = cand[isp]
    jfound[rem[hit]] = j
    keep = np.ones(len(rem), dtype=bool); keep[hit] = False
    rem = rem[keep]
    hist.append(len(rem))
    if j % 50 == 0: print(f'j = {j}: {len(rem):,} classes still without a prime  ({time.time() - t0:.0f} s)', flush=True)
    j += 1
J = jfound[1:]
a_star = int(np.argmax(J) + 1); jmax = int(J.max()); Pq = a_star + jmax * q
wag = (q - 1) * math.log(q) * math.log(q - 1)
lines = [f'q = {q}: least prime found in all {q - 1:,} classes in {time.time() - t0:.0f} s',
         f'P(q) = {Pq:,} (class a = {a_star}, the {jmax + 1}-th member a + {jmax} q); P(q)/q = {Pq / q:.1f}',
         f"Wagstaff's estimate phi(q) ln q ln phi(q) = {wag:,.0f}; P(q) / estimate = {Pq / wag:.3f}; Pomerance's lower order e^gamma phi(q) ln q = {math.exp(0.5772156649) * (q - 1) * math.log(q):,.0f}",
         'classes still without a prime after j members: ' + ', '.join(f'{jj}: {hist[jj - 1]:,}' for jj in (1, 2, 3, 5, 10, 20, 50, 100, 200, 300) if jj - 1 < len(hist))]
top = np.argsort(J)[::-1][:10] + 1
lines.append('the ten latest classes (a: least prime): ' + ', '.join(f'{int(a)}: {int(a) + int(jfound[a]) * q:,}' for a in top))
print('\n'.join(lines))
open(f'probeB_{q}.txt', 'w').write('\n'.join(lines) + '\n')
np.save(f'probeB_{q}_j.npy', jfound.astype(np.int16) if jmax < 32767 else jfound)
