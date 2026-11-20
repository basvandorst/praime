"""Angles of Gaussian primes: all prime ideals of Z[i] with norm <= X (split p = a^2 + b^2: the angles of a + bi and b + ai,
the ramified 1 + i at pi/4, the inert p = 3 mod 4 with p^2 <= X at 0), angles in [0, pi/2).  Number variance in sectors
of a fraction 1/K of [0, pi/2), averaged over 64 positions of the sector grid (approximating Rudnick-Waxman's integral
over theta), divided by the mean N/K.  Rudnick-Waxman Conjecture 1.2: Var ~ (N/K) min(1, 2 log K / log N).
Fakes: lattice points a > 0, b >= 0 with norm n <= X, each kept with chance c / ln n (c chosen so the expected count is N);
and 'odd primitive' fakes: only points with gcd(a, b) = 1 and odd norm.  Prediction written before the run: the primes
follow min(1, 2 log K / log N) roughly, the fakes stay near 1 (minus the thinning, 1 - chance).   usage: rw.py X out.txt [primes]   (primes: skip the fakes, K in steps of 2)"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w'); ONLYP = len(sys.argv) > 3 and sys.argv[3] == 'primes'; STEP = 2 if ONLYP else 4
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
ang = []; latt_n = []; latt_t = []; latt_prim = []
R = int(X ** .5)
for a in range(1, R + 1):
    bmax = int(math.isqrt(X - a * a)); b = np.arange(0, bmax + 1); n = a * a + b * b
    t = np.arctan2(b, a); ok = sv[n]
    ang.append(t[ok])
    if not ONLYP: latt_n.append(n); latt_t.append(t); latt_prim.append((np.gcd(a, b) == 1) & (n % 2 == 1))
inert = [p for p in range(3, int(X ** .5) + 1, 4) if sv[p]]
A = np.concatenate(ang + [np.zeros(len(inert))]); N = len(A); A_split = np.concatenate(ang)
pr('X = %d: N = %d prime ideals (%d split angles incl. 1+i, %d inert at angle 0)' % (X, N, N - len(inert), len(inert)))
if not ONLYP: LN = np.concatenate(latt_n).astype(float); LT = np.concatenate(latt_t); LP = np.concatenate(latt_prim)
def nv(t, K, offs=64):
    w = (math.pi / 2) / K; v = []
    for o in range(offs):
        c = np.bincount(np.floor(((t + o * w / offs) % (math.pi / 2)) / w).astype(int), minlength=K)[:K]
        v.append(c.var())
    return float(np.mean(v)) / (len(t) / K)
rng = np.random.default_rng(20)
if ONLYP:
    pr('    K   log K/log N   primes Var/mean   RW min(1, 2 logK/logN)   without the inert ideals   (inert^2 / N = %.4f)' % (len(inert) ** 2 / N))
    K = 16
    while K <= N // 4:
        bb = math.log(K) / math.log(N); pr('%7d %10.3f %15.3f %19.3f %22.3f' % (K, bb, nv(A, K), min(1, 2 * bb), nv(A_split, K))); K *= STEP
    sys.exit()
c1 = N / float(np.sum(1 / np.log(np.maximum(LN, 2)))); p1 = np.minimum(1, c1 / np.log(np.maximum(LN, 2))); p1[LN < 2] = 0
c2 = N / float(np.sum(LP / np.log(np.maximum(LN, 2)))); p2 = np.where(LP, np.minimum(1, c2 / np.log(np.maximum(LN, 2))), 0.0); p2[LN < 2] = 0
F1 = LT[rng.random(len(LN)) < p1]; F2 = LT[rng.random(len(LN)) < p2]
pr('    K   log K/log N   primes Var/mean   RW min(1, 2 logK/logN)   fakes (all points)   fakes (odd primitive)')
K = 4
while K <= N // 4:
    b = math.log(K) / math.log(N)
    pr('%7d %10.3f %15.3f %19.3f %20.3f %20.3f' % (K, b, nv(A, K), min(1, 2 * b), nv(F1, K), nv(F2, K)))
    K *= 4
