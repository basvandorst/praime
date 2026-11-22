import numpy as np, math, sys, json
import mpmath as mp
X = int(float(sys.argv[1])); SEG = 10**7
mp.mp.dps = 30
chars = {'chi4': (4, [0,1,0,-1]), 'chi3': (3, [0,1,-1])}
cst = {}
for k,(q,chi) in chars.items():
    cst[k] = (np.longdouble(str(mp.dirichlet(1,chi,1)/mp.dirichlet(1,chi))), np.longdouble(str(mp.dirichlet(0,chi,1)/mp.dirichlet(0,chi))))
def small(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, math.isqrt(n)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
sp = small(math.isqrt(X)+1)
# prime powers k>=2 with exact integer arithmetic
pp = []
for p in sp.tolist():
    v = p*p; k = 2
    while v <= X:
        pp.append((v, p, k)); v *= p; k += 1
pp.sort()
ppn = np.array([a for a,_,_ in pp], np.int64); ppp = np.array([b for _,b,_ in pp], np.int64); ppk = np.array([c for _,_,c in pp], np.int64)
def Tfun(y):  # sum_{n odd} y^-n/(n(n+1)), y >= 100
    y = np.asarray(y, np.longdouble); u = 1/y
    return np.arctanh(u) + y/2*np.log1p(-u*u)
CNT = [0]
state = {k: dict(A=np.longdouble(0), B=np.longdouble(0), best=(0.0,None), lastn=None, lasta=None, lastb=None, top_by_dec={}) for k in chars}
def upd(st, r, xv):
    i = int(np.argmax(r))
    if r[i] > st['best'][0]: st['best'] = (float(r[i]), float(xv[i]))
for lo in range(0, X+1, SEG):
    hi = min(lo+SEG, X+1)
    M = np.ones(hi-lo, bool)
    for p in sp.tolist():
        if p*p >= hi: break
        st = max(p*p, -(-lo//p)*p); M[st-lo::p] = False
    if lo == 0: M[:2] = False
    P = np.nonzero(M)[0].astype(np.int64) + lo
    sel = (ppn >= lo) & (ppn < hi)
    n = np.concatenate([P, ppn[sel]]); base = np.concatenate([P, ppp[sel]]); kk = np.concatenate([np.ones(len(P), np.int64), ppk[sel]])
    o = np.argsort(n, kind='stable'); n = n[o]; base = base[o]; kk = kk[o]
    CNT[0] += len(P)
    if len(n) == 0: continue
    lg = np.log(base.astype(np.longdouble)); N = n.astype(np.longdouble)
    for name,(q,chi) in chars.items():
        stt = state[name]; d1, d0 = cst[name]
        cv = np.array(chi, np.longdouble)[base % q] ** kk
        a = cv*lg
        A = stt['A'] + np.cumsum(a/N); B = stt['B'] + np.cumsum(a)
        al = A + d1; be = -B - d0
        # breakpoints: include previous segment's last one so intervals across boundaries are handled
        if stt['lastn'] is not None:
            Nn = np.concatenate([[stt['lastn']], N]); al2 = np.concatenate([[stt['lasta']], al]); be2 = np.concatenate([[stt['lastb']], be])
        else:
            Nn, al2, be2 = N, al, be
        D = al*N + be + Tfun(N)
        m = N >= 100
        if m.any(): upd(stt, np.where(m, np.abs(D)/np.sqrt(N), 0), N)
        # stationary points of (al y + be)/sqrt(y) in (Nn[i], Nn[i+1])
        xs = be2[:-1]/al2[:-1]
        ok = (xs > Nn[:-1]) & (xs < Nn[1:]) & (xs >= 100)
        if ok.any():
            y = xs[ok]; v = np.abs(al2[:-1][ok]*y + be2[:-1][ok] + Tfun(y))/np.sqrt(y); upd(stt, v, y)
        # also y = 100 exactly
        if lo == 0:
            j = np.searchsorted(N, 100, side='right')-1
            y = np.longdouble(100); v = abs(al[j]*y + be[j] + Tfun(y))/10; upd(stt, np.array([v]), np.array([y]))
        stt['A'] = A[-1]; stt['B'] = B[-1]; stt['lastn'] = N[-1]; stt['lasta'] = al[-1]; stt['lastb'] = be[-1]
        # record per-decade top (for comparison)
    if (hi-1) in (10**7, 10**8, 10**9) or hi == X+1:
        print(hi-1, {k: state[k]['best'] for k in chars}, flush=True)
# tail interval from last breakpoint to X
for name in chars:
    stt = state[name]; y = np.longdouble(X); v = abs(stt['lasta']*y + stt['lastb'] + Tfun(y))/np.sqrt(y)
    if v > stt['best'][0]: stt['best'] = (float(v), float(y))
    print(name, 'max ratio on [100,%d]' % X, stt['best'], 'C', float(cst[name][0]-cst[name][1]+np.longdouble(math.log(2))))
print("pi(X) =", CNT[0], "prime powers k>=2:", len(pp))
