"""Shared definitions for the round layouts of prime2.html (N = 500,000).
Angles/radii are recomputed from the page's formulas (prime2.html build/polarLayout/L20/L2P)
and checked against the extracted canvas positions data/<layout>.fx.bin/.fy.bin."""
import numpy as np, os
VIS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(VIS, 'data')
N = 500000
MEXP = np.array([2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,
                 11213,19937,21701,23209,44497,86243,110503,132049,216091])
LAYOUTS = ['sacks','sunflower','polar','wheel','binring','decrings','ternrings','doubling',
           'theodorus','binsquares','bintri','btree']

def sieve(n=N):
    s = np.ones(n+1, bool); s[:2] = False
    for p in range(2, int(n**0.5)+1):
        if s[p]: s[p*p::p] = False
    return s

def load_xy(name):
    fx = np.fromfile(os.path.join(DATA, name+'.fx.bin'), '<f4').astype(np.float64)
    fy = np.fromfile(os.path.join(DATA, name+'.fy.bin'), '<f4').astype(np.float64)
    return fx, fy

def bl(n):  # bit length
    return np.floor(np.log2(n)).astype(np.int64) + 1

def layout(name):
    """Return dict with n (1..N), theta (layout angle in [0,2pi) or offset*2pi), ring/radius,
    model canvas x,y (from the page's formula), plus 't' (fraction along ring) where relevant."""
    n = np.arange(1, N+1, dtype=np.float64)
    out = {}
    if name == 'sacks':
        R = int(np.ceil(np.sqrt(N))) + 2
        r = np.sqrt(n); a = 2*np.pi*r
        x, y = R + r*np.cos(a), R - r*np.sin(a)
        out.update(theta=np.mod(a, 2*np.pi), rad=r)
    elif name in ('sunflower', 'polar', 'wheel', 'binring', 'doubling'):
        C = 702.0
        if name == 'sunflower':
            ga = np.pi*(3-np.sqrt(5)); c = 696/np.sqrt(N); r = c*np.sqrt(n); a = n*ga
        elif name == 'polar':
            c = 696/N; r = c*n; a = n.copy()
        elif name == 'wheel':
            w = 708; rings = int(np.ceil(N/w)); c = 666/rings
            ni = np.arange(1, N+1)
            r = 30 + c*((ni-1)//w); a = np.pi/2 - 2*np.pi*((ni-1) % w)/w
            out['ring'] = (ni-1)//w
        elif name == 'binring':
            B = int(np.floor(np.log2(N))) + 1; c = 676/max(1, B-1)
            b = bl(n); h = 2.0**(b-1); t = (n-h+0.5)/h
            r = 20 + c*(b-1); a = np.pi/2 - 2*np.pi*t
            out.update(ring=b, t=t)
        elif name == 'doubling':
            c = 681/np.log2(N); L = np.log2(n); r = 15 + c*L; a = np.pi/2 - 2*np.pi*L
            out['t'] = np.mod(L, 1.0)
        x, y = C + r*np.cos(a), C - r*np.sin(a)
        out.update(theta=np.mod(a, 2*np.pi), rad=r)
    elif name in ('decrings', 'ternrings', 'theodorus'):
        if name == 'theodorus':
            th = np.zeros(N+1)
            th[2:] = np.cumsum(np.arctan(1/np.sqrt(np.arange(1, N, dtype=np.float64))))
            c = 696/np.sqrt(N); r = c*np.sqrt(n); a = th[1:]
        else:
            B = 10 if name == 'decrings' else 3
            Dm = 1
            while B**Dm <= N: Dm += 1
            c = 676/max(1, Dm-1)
            ni = np.arange(1, N+1)
            d = np.ones(N, np.int64)
            for k in range(1, 30):
                d[ni >= B**k] = k+1
            h = np.power(float(B), d-1); cnt = (B-1)*h; t = (n-h+0.5)/cnt
            r = 20 + c*(d-1); a = np.pi/2 - 2*np.pi*t
            out.update(ring=d, t=t)
        x, y = 2 + 700 + r*np.cos(a), 2 + 700 - r*np.sin(a)
        out.update(theta=np.mod(a, 2*np.pi), rad=r)
    elif name in ('binsquares', 'bintri', 'btree'):
        Bm = int(bl(np.array([N]))[0])
        b = bl(n); h = 2.0**(b-1); t = (n-h+0.5)/h
        out.update(ring=b, t=t, rad=b.astype(float), theta=2*np.pi*t)
        if name == 'btree':
            x = 2 + t*1400; y = 2 + (b-1)*1400/max(1, Bm-1)
        elif name == 'bintri':
            half = 700*b/Bm; x = 2 + 700 + (t-0.5)*2*half; y = 2 + (b-1)*1400/max(1, Bm-1)
        else:
            c = 676/max(1, Bm-1); r = 20 + c*(b-1); tt = np.mod(t+0.125, 1.0); s = 8*r*tt
            x = np.empty(N); y = np.empty(N)
            m1 = s < 2*r; m2 = (s >= 2*r) & (s < 4*r); m3 = (s >= 4*r) & (s < 6*r); m4 = s >= 6*r
            x[m1] = -r[m1]+s[m1]; y[m1] = -r[m1]
            x[m2] = r[m2]; y[m2] = -r[m2]+(s[m2]-2*r[m2])
            x[m3] = r[m3]-(s[m3]-4*r[m3]); y[m3] = r[m3]
            x[m4] = -r[m4]; y[m4] = r[m4]-(s[m4]-6*r[m4])
            x = 2 + 700 + x; y = 2 + 700 + y
    out['x'] = np.concatenate([[np.nan], x]); out['y'] = np.concatenate([[np.nan], y])
    for k in list(out):
        if k not in ('x', 'y'):
            out[k] = np.concatenate([[np.nan], np.asarray(out[k], dtype=np.float64)])
    return out

# ---------- null model: pseudo-exponents ----------
def null_weights(isp):
    P = np.nonzero(isp)[0]
    g = np.exp(np.euler_gamma)
    a = np.where(P % 4 == 1, 6.0, 2.0)
    w = np.minimum(1.0, g*np.log2(a*P)/P)
    return P, w

def sample_null(P, w, nsets, m=31, seed=1):
    """Exact sampling of independent Bernoulli(w_p) picks conditioned on exactly m picks
    (backward DP of the Poisson-binomial tail, then a forward sequential draw)."""
    rng = np.random.default_rng(seed)
    L = len(P)
    # Q[i, j] = P(exactly j picks among items i..L-1), j = 0..m
    Q = np.zeros((L+1, m+1)); Q[L, 0] = 1.0
    for i in range(L-1, -1, -1):
        Q[i] = (1-w[i])*Q[i+1]
        Q[i, 1:] += w[i]*Q[i+1, :-1]
        # no renormalisation needed: values stay well above underflow for m = 31
    need = np.full(nsets, m)
    picks = np.zeros((nsets, m), np.int64); cnt = np.zeros(nsets, np.int64)
    U = rng.random((nsets,))
    for i in range(L):
        act = need > 0
        if not act.any(): break
        j = need
        num = w[i]*Q[i+1, np.maximum(j-1, 0)]
        den = Q[i, j]
        pr = np.where(act, num/np.where(den > 0, den, 1), 0.0)
        u = rng.random(nsets)
        take = (u < pr) & act
        idx = np.nonzero(take)[0]
        picks[idx, cnt[idx]] = P[i]; cnt[idx] += 1; need[idx] -= 1
    assert (cnt == m).all()
    return picks
