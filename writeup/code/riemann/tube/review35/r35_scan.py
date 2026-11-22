"""Own scan: max over y in [Y0, X] of |G(y) - ln y - 1 - B0 + T_e(y)|/sqrt(y), chi = (./5).
G piecewise: on [n_i, n_{i+1}] G = alpha x - B_i, alpha = A_i + d1.  f = alpha x + beta - ln x (+T_e), beta = -B_i - 1 - B0.
Interior extrema of f/sqrt(x): h(x) = alpha x - beta + ln x - 2 = 0 (T_e neglected in h, included in values); solved by bisection
on monotone pieces (split at x = -1/alpha). Accumulation in long double."""
import sys, math, json
import numpy as np
X = int(float(sys.argv[1])); Y0 = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
d1 = 0.8276794767155048879910469896773736467539; B0 = -0.02202465783872691842463476000992357551451
C = 0.1565569539942864969984496282291206541929
LD = np.longdouble
def Te(x):
    u = 1.0 / x
    return -0.5 * np.log1p(-u * u) - np.arctanh(u) / u + 1.0 if np.ndim(x) else (-0.5*math.log1p(-u*u) - math.atanh(u)/u + 1.0 if x > 1 else 1 - math.log(2))
r = math.isqrt(X) + 1
bs = np.ones(r + 1, bool); bs[:2] = False
for i in range(2, math.isqrt(r) + 1):
    if bs[i]: bs[i*i::i] = False
base = np.nonzero(bs)[0]
pw = []
for p in base.tolist():
    v = p * p; k = 2
    while v <= X:
        pw.append((v, p, k)); v *= p; k += 1
pw.sort(); pwn = np.array([t[0] for t in pw], np.int64)
pwa = np.array([math.log(t[1]) * ([0, 1, -1, -1, 1][t[1] % 5] ** t[2]) for t in pw])
chi = np.array([0, 1, -1, -1, 1], float)
bi = (0.0, None); A = LD(0); B = LD(0); best = (0.0, None, ''); seg = 10**7
top = []  # top few local maxima
prev_last = None  # (x, alpha, beta) of last breakpoint of previous segment
def evalg(x, al, be):
    return np.abs(al * x + be - np.log(x) + Te(x)) / np.sqrt(x)
for s0 in range(0, X + 1, seg):
    s1 = min(s0 + seg, X + 1)
    m = np.ones(s1 - s0, bool)
    if s0 == 0: m[:2] = False
    for p in base.tolist():
        if p * p >= s1: break
        st = max(p * p, ((s0 + p - 1) // p) * p)
        m[st - s0::p] = False
    P = np.nonzero(m)[0].astype(np.int64) + s0
    if len(P) == 0 and not ((pwn >= s0) & (pwn < s1)).any(): continue
    a = np.log(P.astype(float)) * chi[P % 5]
    sel = (pwn >= s0) & (pwn < s1)
    nn = np.concatenate([P, pwn[sel]]); aa = np.concatenate([a, pwa[sel]])
    o = np.argsort(nn, kind='stable'); nn = nn[o]; aa = aa[o]
    cA = A + np.cumsum(aa.astype(LD) / nn.astype(LD)); cB = B + np.cumsum(aa.astype(LD))
    A = cA[-1]; B = cB[-1]
    x = nn.astype(LD)
    al = cA + LD(d1); be = -cB - LD(1 + B0)
    # breakpoint values
    v = evalg(x.astype(float), al, be)  # mixed: al,be long double
    v = np.asarray(v, float)
    msk = x >= Y0
    if msk.any():
        i = int(np.argmax(np.where(msk, v, -1)))
        if v[i] > best[0]: best = (float(v[i]), float(x[i]), 'bp')
    # intervals [x_i, x_{i+1}] ; include the interval bridging from previous segment
    xl = x; xr = np.concatenate([x[1:], [LD(np.inf)]])
    if prev_last is not None:
        xl = np.concatenate([[prev_last[0]], xl]); xr = np.concatenate([[x[0]], xr])
        al2 = np.concatenate([[prev_last[1]], al]); be2 = np.concatenate([[prev_last[2]], be])
    else:
        al2, be2 = al, be
    prev_last = (x[-1], al[-1], be[-1])
    xl = xl[:-1]; xr = xr[:-1]; al2 = al2[:-1]; be2 = be2[:-1]
    xl = np.maximum(xl, LD(Y0)); ok = xr > xl
    xl, xr, al2, be2 = xl[ok], xr[ok], al2[ok], be2[ok]
    # boundary point Y0 itself
    for pieces in range(2):
        # split at xc = -1/alpha
        xc = np.where(al2 < 0, -1 / al2, LD(-1))
        if pieces == 0: lo, hi = xl, np.where((xc > xl) & (xc < xr), xc, xr)
        else:
            inn = (xc > xl) & (xc < xr)
            lo, hi = xc[inn], xr[inn]; alp, bep = al2[inn], be2[inn]
        if pieces == 0: alp, bep = al2, be2
        h = lambda t: alp * t - bep + np.log(t) - 2  # long double ok
        hl, hh = h(lo), h(hi)
        cand = np.sign(hl) != np.sign(hh)
        lo_, hi_, a_, b_ = lo[cand], hi[cand], alp[cand], bep[cand]
        hl_ = hl[cand]
        for _ in range(80):
            mid = (lo_ + hi_) / 2; hm = a_ * mid - b_ + np.log(mid) - 2
            same = np.sign(hm) == np.sign(hl_)
            lo_ = np.where(same, mid, lo_); hi_ = np.where(same, hi_, mid); hl_ = np.where(same, hm, hl_)
        xs = (lo_ + hi_) / 2
        if len(xs):
            vs = np.asarray(evalg(xs.astype(float), a_, b_), float)
            j = int(np.argmax(vs))
            if vs[j] > bi[0]: bi = (float(vs[j]), float(xs[j]))
            if vs[j] > best[0]: best = (float(vs[j]), float(xs[j]), 'interior')
    # also boundary Y0 and X
    for xb in (Y0, float(X)):
        if s0 <= xb < s1 or (xb == X and s1 == X + 1):
            k = np.searchsorted(nn, int(xb), side='right') - 1
            if k >= 0:
                vb = float(evalg(xb, al[k], be[k]))
                if vb > best[0]: best = (vb, xb, 'boundary')
    print(s1, best, flush=True) if s0 % (10 * seg) == 0 else None
print(json.dumps({'X': X, 'Y0': Y0, 'max': best[0], 'at': best[1], 'kind': best[2], 'ratio_C': best[0] / C, 'best_interior': bi}))
