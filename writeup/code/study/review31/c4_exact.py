"""Review 31: exact I(X, d) = int_1^X (psi((1+d)x) - psi(x) - d x)^2 dx (piecewise-exact), own segmented sieve.
usage: c4_exact.py delta Xmax  -> R at the r96 X grid points <= Xmax, saved to c4_<delta>.json"""
import sys, math, json, time
import numpy as np

d = float(sys.argv[1]); XMAX = float(sys.argv[2]); q = 1.0 + d
Xs_all = np.exp(np.linspace(math.log(1e7), math.log(1e10), 200))
Xs = Xs_all[Xs_all <= XMAX] if len(sys.argv) < 4 else np.array([float(sys.argv[3])])
TOP = int(q * XMAX) + 10
base = None

def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]

BP = small_primes(int(math.isqrt(TOP)) + 2)

def primes_in(lo, hi):
    """primes in [lo, hi)"""
    if hi <= lo: return np.zeros(0, np.int64)
    s = np.ones(hi - lo, bool)
    for p in BP.tolist():
        if p * p >= hi: break
        st = max(p * p, ((lo + p - 1) // p) * p)
        s[st - lo::p] = False
    if lo < 2: s[:2 - lo] = False
    return np.nonzero(s)[0].astype(np.int64) + lo

# prime powers k >= 2
pp = []
for p in BP.tolist():
    v = p * p
    while v <= TOP:
        pp.append((v, math.log(p))); v *= p
pp.sort(); PPV = np.array([a for a, _ in pp], np.int64); PPL = np.array([b for _, b in pp])

def jumps(lo, hi):
    """prime powers n in [lo, hi) and their Lambda"""
    P = primes_in(lo, hi)
    m = (PPV >= lo) & (PPV < hi)
    v = np.concatenate([P, PPV[m]]); l = np.concatenate([np.log(P.astype(np.float64)), PPL[m]])
    o = np.argsort(v, kind='stable'); return v[o], l[o]

t0 = time.time()
L = 10_000_000
s = 1; A = 0.0; Itot = []; out = {}
Xi = 0
while s < XMAX:
    e = min(s + L, int(XMAX) + 2)
    n1, l1 = jumps(s, e)                                   # psi(x) jumps at x = n: A -= Lambda
    M0 = math.ceil(q * s); M1 = math.ceil(q * e)
    n2, l2 = jumps(M0, M1)                                 # psi(qx) jumps at x = m/q: A += Lambda
    x2 = np.clip(n2 / q, s, e)
    xm = Xs[(Xs >= s) & (Xs < e)]
    pos = np.concatenate([n1.astype(np.float64), x2, xm])
    dl = np.concatenate([-l1, l2, np.zeros(len(xm))])
    tag = np.concatenate([np.zeros(len(n1) + len(x2), np.int8), np.ones(len(xm), np.int8)])
    o = np.argsort(pos, kind='stable'); pos = pos[o]; dl = dl[o]; tag = tag[o]
    # intervals: [s, pos0), [pos0, pos1), ..., [pos_last, e)
    left = np.concatenate([[float(s)], pos]); right = np.concatenate([pos, [float(e)]])
    Aval = A + np.concatenate([[0.0], np.cumsum(dl)])
    y = Aval - d * left; D = right - left
    piece = y * y * D - y * d * D * D + d * d * D ** 3 / 3.0
    cs = np.cumsum(piece)
    base_I = Itot[-1] if Itot else 0.0
    # I at X marks: integral up to pos of the mark = cs before the interval starting at that mark
    idx = np.nonzero(tag)[0]                               # mark k is pos[k]; integral up to it = sum piece[0..k]
    for k in idx.tolist():
        X = pos[k]; I = base_I + float(cs[k])
        out[repr(float(X))] = I / (0.5 * d * X * X)
    Itot.append(base_I + float(cs[-1]))
    A = float(Aval[-1]); s = e
    print('seg to %.3e  A-?  %.1fs' % (e, time.time() - t0), flush=True) if (e // L) % 20 == 0 else None
R = np.array([out[repr(float(X))] for X in Xs])
json.dump({'delta': d, 'X': Xs.tolist(), 'R': R.tolist(), 'seconds': time.time() - t0}, open('/home/bas/scratch/study-review31/c4_%g_%g.json' % (d, XMAX), 'w'))
print('done', len(Xs), 'avg R %.5f' % R.mean(), time.time() - t0)
