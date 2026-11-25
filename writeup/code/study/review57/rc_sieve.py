"""Own recount (review of round 155).  Segmented odd-only numpy sieve on (A, B].
Per cell (X[i-1], X[i]] of the grid X_k = floor(exp(k/2048)) it accumulates, for primes p in the cell:
  th      sum ln p
  N2,P2   twin p (p+2 prime): count, sum ln p * ln(p+2)
  N4,P4   cousin p (p+4 prime)
  NSG,SG  Sophie Germain p (2p+1 prime): count, sum ln p * ln(2p+1)
  thc     sum ln p by p mod 420   (gives classes mod 3,4,5,7,105,...)
  twc     sum ln p ln(p+2) by p mod 420 (twins)
  cuc     sum ln p ln(p+4) by p mod 420 (cousins)
Definitions are from the write-up text (twinzeros.md, section 3), not from the authors' code.
usage: rc_sieve.py hold.npz out.npz i0 i1      (cells i0+1 .. i1, i.e. p in (X[i0], X[i1]])
"""
import sys, time, math
import numpy as np

hold, out, i0, i1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
X = np.load(hold)['X'].astype(np.int64)
A, B = int(X[i0]), int(X[i1])
print('range (A,B] =', A, B, 'length', B - A, 'cells', i1 - i0, flush=True)

def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]

R1 = small_primes(math.isqrt(B + 6) + 1)          # for p, p+2, p+4, p+6
R2 = small_primes(math.isqrt(2 * B + 1) + 1)      # for 2p+1
R1 = R1[R1 > 2]; R2 = R2[R2 > 2]
M = 420
nX = len(X)
th = np.zeros(nX); N2 = np.zeros(nX); P2 = np.zeros(nX); N4 = np.zeros(nX); P4 = np.zeros(nX)
NSG = np.zeros(nX); SG = np.zeros(nX); NT026 = np.zeros(nX); T026 = np.zeros(nX); NT046 = np.zeros(nX); T046 = np.zeros(nX)
thc = np.zeros(nX * M); twc = np.zeros(nX * M); cuc = np.zeros(nX * M)

L = 1 << 23                     # odd entries per segment
s0 = A + 1 if (A + 1) % 2 == 1 else A + 2   # first odd > A
t0 = time.time()
seg = s0
nseg = 0
while seg <= B:
    n = min(L, (B - seg) // 2 + 1)          # odd numbers seg, seg+2, ..., in (A,B]
    ne = n + 3                                # extra entries for p+2, p+4, p+6
    a = np.ones(ne, dtype=np.uint8)           # a[j] = 1 iff seg+2j has no odd prime factor r <= sqrt
    b = np.ones(n, dtype=np.uint8)            # b[j] = 1 iff 2(seg+2j)+1 has no prime factor r <= sqrt(2B+1)
    for r in R1.tolist():
        inv2 = (r + 1) // 2
        j = ((-seg) * inv2) % r
        a[j::r] = 0
    for r in R2.tolist():
        inv2 = (r + 1) // 2
        t = (r - 1) // 2
        j = ((t - seg) * inv2) % r
        b[j::r] = 0
    idx = np.nonzero(a[:n])[0]
    p = seg + 2 * idx.astype(np.int64)
    lp = np.log(p.astype(np.float64))
    cell = np.searchsorted(X, p, side='left')          # X[cell-1] < p <= X[cell]
    cls = p % M
    th += np.bincount(cell, weights=lp, minlength=nX)
    thc += np.bincount(cell * M + cls, weights=lp, minlength=nX * M)
    m2 = a[idx + 1].astype(bool)
    w2 = lp[m2] * np.log((p[m2] + 2).astype(np.float64))
    N2 += np.bincount(cell[m2], minlength=nX); P2 += np.bincount(cell[m2], weights=w2, minlength=nX)
    twc += np.bincount(cell[m2] * M + cls[m2], weights=w2, minlength=nX * M)
    m4 = a[idx + 2].astype(bool)
    w4 = lp[m4] * np.log((p[m4] + 4).astype(np.float64))
    N4 += np.bincount(cell[m4], minlength=nX); P4 += np.bincount(cell[m4], weights=w4, minlength=nX)
    cuc += np.bincount(cell[m4] * M + cls[m4], weights=w4, minlength=nX * M)
    ms = b[idx].astype(bool)
    ws = lp[ms] * np.log((2 * p[ms] + 1).astype(np.float64))
    NSG += np.bincount(cell[ms], minlength=nX); SG += np.bincount(cell[ms], weights=ws, minlength=nX)
    m6 = a[idx + 3].astype(bool)
    l2 = np.log((p + 2).astype(np.float64)); l4 = np.log((p + 4).astype(np.float64)); l6 = np.log((p + 6).astype(np.float64))
    t1 = m2 & m6; t2 = m4 & m6
    NT026 += np.bincount(cell[t1], minlength=nX); T026 += np.bincount(cell[t1], weights=lp[t1] * l2[t1] * l6[t1], minlength=nX)
    NT046 += np.bincount(cell[t2], minlength=nX); T046 += np.bincount(cell[t2], weights=lp[t2] * l4[t2] * l6[t2], minlength=nX)
    seg += 2 * n
    nseg += 1
    if nseg % 20 == 0:
        print('seg', nseg, 'at', seg, '%.1fs' % (time.time() - t0), flush=True)
np.savez_compressed(out, X=X, i0=i0, i1=i1, th=th, N2=N2, P2=P2, N4=N4, P4=P4, NSG=NSG, SG=SG, NT026=NT026, T026=T026, NT046=NT046, T046=T046,
         thc=thc.reshape(nX, M), twc=twc.reshape(nX, M), cuc=cuc.reshape(nX, M))
print('done %.1fs' % (time.time() - t0), 'primes', 'twins', int(N2.sum()), 'SG', int(NSG.sum()))
