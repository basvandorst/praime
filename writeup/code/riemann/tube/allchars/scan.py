"""max_{100<=y<=X} |G_chi(y) - c(y)|/sqrt(y) for all characters in consts.json plus zeta.
G piecewise linear in y between prime powers (for zeta / even chi the target c(y) has a ln y or y ln y term).
On a piece G - c = alpha*y + beta with constant complex alpha, beta (odd chi): |.|^2/y = |alpha|^2 y + 2Re(alpha conj beta) + |beta|^2/y
is convex in y, so the max over a piece is at an endpoint. Interior points are sampled for y < 1e6 as a check for all chars.
usage: scan.py X"""
import sys, json, math
import numpy as np

X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 8
cons = json.load(open('consts.json'))
from chars import chars
EG = 0.5772156649015329
LN2PI = math.log(2 * math.pi)


def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


base = small_primes(int(math.isqrt(X)) + 2)
pw = []   # (v, p, k) k >= 2
for p in base.tolist():
    v, k = p * p, 2
    while v <= X:
        pw.append((v, p, k)); v *= p; k += 1
pw.sort()
pwv = np.array([t[0] for t in pw], np.int64); pwp = np.array([t[1] for t in pw], np.int64); pwk = np.array([t[2] for t in pw], np.int64)

C = []
for name, q, v, a, real in chars():
    c = cons[name]
    C.append(dict(name=name, q=q, chi=np.array(v, complex), a=a, dl1=complex(*c['dl1']), c0=complex(*c['c0']), kind='odd' if a else 'even'))
C.append(dict(name='zeta', q=1, chi=np.array([1.0 + 0j]), a=0, dl1=0j, c0=complex(LN2PI), kind='zeta'))
for c in C:
    c.update(A=0j, B=0j, best=0.0, at=None, phase=None, ibest=0.0, iat=None, dec=np.zeros(10))


def target(c, x):
    if c['kind'] == 'even':
        return np.log(x) + 1 + c['c0']
    return c['c0'] + 0 * x


def Gval(c, x, A, B):
    if c['kind'] == 'zeta':
        return x * A - B - x * np.log(x) + (1 + EG) * x
    return x * (A + c['dl1']) - B


seg = 10 ** 7
for s0 in range(0, X + 1, seg):
    s1 = min(s0 + seg, X + 1)          # [s0, s1)
    lo = max(s0, 2)
    mask = np.ones(s1 - lo, bool)
    for p in base.tolist():
        if p * p >= s1:
            break
        st = max(p * p, ((lo + p - 1) // p) * p)
        mask[st - lo::p] = False
    P = np.nonzero(mask)[0].astype(np.int64) + lo
    sel = (pwv >= s0) & (pwv < s1)
    nn = np.concatenate([P, pwv[sel]]); pp = np.concatenate([P, pwp[sel]]); kk = np.concatenate([np.ones(len(P), np.int64), pwk[sel]])
    w = np.log(pp.astype(float))
    for extra in (100, X):                # evaluation points with weight 0
        if s0 <= extra < s1:
            nn = np.append(nn, extra); pp = np.append(pp, 1); kk = np.append(kk, 1); w = np.append(w, 0.0)
    o = np.argsort(nn, kind='stable'); nn, pp, kk, w = nn[o], pp[o], kk[o], w[o]
    x = nn.astype(float)
    nxt = np.append(x[1:], float(s1))
    for c in C:
        chv = c['chi'][pp % c['q']] ** kk
        an = w * chv
        cA = c['A'] + np.cumsum(an / x); cB = c['B'] + np.cumsum(an)
        G = Gval(c, x, cA, cB)
        r = np.abs(G - target(c, x)) / np.sqrt(x)
        m = (x >= 100) & (x <= X)
        rr = np.where(m, r, 0.0)
        i = int(np.argmax(rr))
        if rr[i] > c['best']:
            c['best'] = float(rr[i]); c['at'] = float(x[i]); c['phase'] = float(np.angle(G[i] - target(c, x[i])))
        d = np.floor(np.log10(x[m])).astype(int)
        np.maximum.at(c['dec'], d, r[m])
        # interior samples for y < 1e6
        mi = m & (x < 1e6) & (nxt <= X)
        if mi.any():
            for t in np.linspace(0.05, 0.95, 19):
                y = x[mi] + t * (nxt[mi] - x[mi])
                Gy = Gval(c, y, cA[mi], cB[mi])
                ry = np.abs(Gy - target(c, y)) / np.sqrt(y)
                j = int(np.argmax(ry))
                if ry[j] > c['ibest']:
                    c['ibest'] = float(ry[j]); c['iat'] = float(y[j])
        c['A'] = complex(cA[-1]); c['B'] = complex(cB[-1])
    print('done', s1, flush=True)

res = {}
for c in C:
    C_ = cons[c['name']]['Cb'] if c['name'] != 'zeta' else 2 + EG - math.log(4 * math.pi)
    res[c['name']] = dict(q=c['q'], kind=c['kind'], C=C_, max=c['best'], at=c['at'], frac=c['best'] / C_, phase=c['phase'],
                          interior_max_below_1e6=c['ibest'], interior_at=c['iat'], decade_max=c['dec'].tolist())
    print('%-16s q=%2d %-5s C=%.6f max=%.6f at %.0f frac=%.3f  interior(<1e6) %.6f at %.1f' % (
        c['name'], c['q'], c['kind'], C_, c['best'], c['at'], c['best'] / C_, c['ibest'], c['iat']))
json.dump(res, open('scan_%d.json' % X, 'w'), indent=1)
