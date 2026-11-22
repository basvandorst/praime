"""Independent recomputation of R(X, delta) = int_1^X (psi((1+d)x) - psi(x) - d x)^2 dx / (d X^2 / 2).
Own odd-only segmented sieve; psi via sorted prime table + cumulative logs; prime powers added.
mode 'mid NPT': midpoint rule with NPT points; mode 'exact': exact piecewise integration (step function)."""
import sys, math, time
import numpy as np

X = int(float(sys.argv[1])); mode = sys.argv[2]
DELTAS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
YMAX = int(X * 1.1) + 10
t0 = time.time()

def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

bp = base_primes(math.isqrt(YMAX) + 1)
chunks = [np.array([2], np.uint32)]
SEG = 1 << 24  # odd numbers per segment
lo = 1  # segment covers odd numbers lo, lo+2, ..., lo + 2(SEG-1)
while lo <= YMAX:
    hi = min(lo + 2 * SEG, YMAX + 1)  # exclusive
    m = (hi - lo + 1) // 2
    seg = np.ones(m, bool)
    for p in bp[1:]:
        p = int(p)
        if p * p >= hi: break
        st = max(p * p, ((lo + p - 1) // p) * p)
        if st % 2 == 0: st += p
        seg[(st - lo) // 2::p] = False
    if lo == 1: seg[0] = False
    chunks.append((np.nonzero(seg)[0] * 2 + lo).astype(np.uint32))
    lo = hi if hi % 2 == 1 else hi + 1
P = np.concatenate(chunks); del chunks
P = P[P <= YMAX]
CUM = np.empty(len(P) + 1); CUM[0] = 0.0; base = 0.0; B = 1 << 22
for i in range(0, len(P), B):
    L = np.log(P[i:i+B].astype(np.float64))
    c = np.cumsum(L); CUM[i+1:i+1+len(c)] = base + c; base += math.fsum(L.tolist())
del L
pw = []
for p in bp.tolist():
    v = p * p
    while v <= YMAX:
        pw.append((v, math.log(p))); v *= p
pw.sort(); PWV = np.array([a for a, _ in pw], float); PWC = np.concatenate([[0.0], np.cumsum([b for _, b in pw])])
print('sieve: %d primes to %d, %.1f s' % (len(P), YMAX, time.time() - t0), flush=True)

def psi(q):
    fq = np.floor(q)
    k = np.searchsorted(P, fq, side='right')
    j = np.searchsorted(PWV, fq, side='right')
    return CUM[k] + PWC[j]

print('psi(1e8)=%.4f' % psi(np.array([1e8]))[0], ' psi(X)=%.4f' % psi(np.array([float(X)]))[0], flush=True)

res = {}
if mode == 'mid':
    NPT = int(float(sys.argv[3]))
    h = (X - 1) / NPT
    acc = {d: 0.0 for d in DELTAS}
    CH = 2_000_000
    for i0 in range(0, NPT, CH):
        x = (np.arange(i0, min(NPT, i0 + CH)) + 0.5) * h + 1.0
        p0 = psi(x)
        for d in DELTAS:
            f = psi(x * (1 + d)) - p0 - d * x
            acc[d] += float(np.dot(f, f)) * h
    for d in DELTAS: res[d] = acc[d] / (0.5 * d * X * X)
else:
    # exact: on each interval between consecutive breakpoints, A = psi((1+d)x) - psi(x) is constant;
    # int_a^b (A - d x)^2 dx = (v - u)(u^2 + u v + v^2) / (3 d) with u = d a - A, v = d b - A.
    dl = [float(a) for a in sys.argv[3].split(',')] if len(sys.argv) > 3 else DELTAS
    for d in dl:
        acc = 0.0
        CHX = 10_000_000
        for a0 in range(1, X, CHX):
            b0 = min(X, a0 + CHX)
            def rng(lo_, hi_):
                r = [P[np.searchsorted(P, lo_, 'right'):np.searchsorted(P, hi_, 'left')].astype(np.float64)]
                r.append(PWV[(PWV > lo_) & (PWV < hi_)])
                return np.concatenate(r)
            br1 = rng(a0, b0); br2 = rng(a0 * (1 + d), b0 * (1 + d)) / (1 + d)
            br = np.concatenate([[float(a0)], br1, br2, [float(b0)]]); br.sort()
            a = br[:-1]; b = br[1:]; mid = 0.5 * (a + b)
            A = psi(mid * (1 + d)) - psi(mid)
            u = d * a - A; v = d * b - A
            acc += float(np.sum((v - u) * (u * u + u * v + v * v))) / (3 * d)
        res[d] = acc / (0.5 * d * X * X)
        print('  exact delta %.0e R %.5f  (%.0f s)' % (d, res[d], time.time() - t0), flush=True)
Bc = 1 - 0.5772156649015329 - math.log(2 * math.pi)
for d in sorted(res, reverse=True):
    print('%s X=%.0e delta %.0e  R %.5f   R-(ln(1/d)+B) %+.4f' % (mode, X, d, res[d], res[d] - math.log(1 / d) - Bc))
print('total %.0f s' % (time.time() - t0))
