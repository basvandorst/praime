# Explicit zero exclusion from a one-sided bound H(y) <= B(y), theta = 1/2 (vectorised, log space).
import math, sys, numpy as np, mpmath as mp
from scipy.special import gammainc, gammaln
W = 0.05; env = np.load('env1e8.npy'); B1 = math.log(1e8)
CRH = 2 + 0.57721566490153286061 - math.log(4*math.pi)   # 0.0461914
GAM0 = 14.0
def Nerr(T): return 0.112*math.log(T) + 0.278*math.log(math.log(T)) + 2.51   # Trudgian 2014
def Nmain(T): return T/(2*math.pi)*math.log(T/(2*math.pi*math.e)) + 7/8
def Nup(T): return 0.0 if T < GAM0 else Nmain(T) + Nerr(T)
S1 = 2*float(mp.quad(lambda t: Nup(float(t))*2*t/(0.25+t*t)**2, [GAM0, 100, 1e4, 1e6, mp.inf]))
def nmax(g0, D):
    lo = g0 - D
    if lo <= GAM0: return max(1, math.floor(Nup(g0 + D)))
    return max(1, math.floor(Nmain(g0+D) + Nerr(g0+D) - Nmain(lo) + Nerr(lo)))
def logQ(nu, x):   # log of regularized upper gamma Q(nu,x), rigorous upper bound
    with np.errstate(all='ignore'):
        bd = (nu-1)*np.log(x) - x - gammaln(nu) - np.log1p(-(nu-1)/x)
    return np.where(x > nu - 1 + 1e-9, np.minimum(bd, 0.0), 0.0)
def logP(nu, x):   # log of regularized lower gamma P(nu,x), upper bound
    with np.errstate(all='ignore'):
        bd = nu*np.log(x) - x - gammaln(nu+1) - np.log1p(-x/(nu+1))
    return np.where(x < nu + 1 - 1e-9, np.minimum(bd, 0.0), 0.0)
def lse(*xs): 
    M = np.maximum.reduce(xs); return M + np.log(sum(np.exp(x - M) for x in xs))
ub = np.arange(len(env))*W; ue = np.minimum(ub + W, B1)
_cache = {}
def upper(lam, mode, b, NU):
    key = (lam, mode, b, len(NU))
    if key in _cache: return _cache[key]
    if len(_cache) > 40: _cache.clear()
    sig = 1 + lam; R = sig - 0.5; t = []
    if mode == 'full':
        Bp = np.maximum(env, 0); J = np.zeros(len(NU))
        for i0 in range(0, len(NU), 2000):
            nu = NU[i0:i0+2000, None]
            J[i0:i0+2000] = ((gammainc(nu, R*ue[None, :]) - gammainc(nu, R*ub[None, :]))*Bp[None, :]).sum(1)
        with np.errstate(divide='ignore'):
            lJ = np.where(J > 1e-280, np.log(np.maximum(J, 1e-300)), math.log(Bp.max()) + logP(NU, R*b))
        t.append(math.log(2) + lJ - NU*math.log(R))
    elif mode == 'hyp':      # H <= C sqrt(y) on [1, e^b]
        t.append(math.log(2*CRH) + logP(NU, R*b) - NU*math.log(R))
    else:                    # window [e^a, e^b] with C0, trivial below
        a, C0 = mode
        t.append(math.log(2*C0) - NU*math.log(R))
        t.append(math.log(2*S1) + logP(NU, lam*a) - NU*math.log(lam))
    t.append(math.log(2*S1) + logQ(NU, lam*b) - NU*math.log(lam))
    t.append(math.log(S1) - NU*0.5*math.log(lam**2 + GAM0**2))
    _cache[key] = lse(*t); return _cache[key]
LAMS = [0.02,0.05,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,4,6,10,15,25,40,60,100,150,250,400,700,1000]
NUMAX = 400000
def excluded(beta, g0, mode, rig, b=B1):
    rho = complex(beta, g0); r0 = 1/abs(rho*(1-rho))
    for lam in LAMS:
        sig = 1 + lam; d0 = sig - beta; R = sig - 0.5
        nmx = min(NUMAX if mode != "full" else 40000, int(3*lam*b) + 50); NU = np.arange(1, nmx + 1, dtype=float)
        U = upper(lam, mode, b, NU)
        if not rig:
            main = math.log(r0) - NU*math.log(d0)
            far = math.log(CRH + 3*r0) - NU*math.log(R)
            with np.errstate(invalid='ignore', divide='ignore'):
                low = main + np.log1p(-np.exp(np.minimum(far - main, 0)))
            ok = (far < main) & (low > U)
            if ok.any(): return (lam, int(NU[ok.argmax()]))
        else:
            for Df in [1.02, 1.05, 1.1, 1.2, 1.4, 1.7, 2.0, 3.0]:
                D = d0*Df; n = nmax(g0, D); lrmin = math.log((g0+D)**2/(1+(g0+D)**2)**2)
                for m in range(0, nmx - n, max(1, nmx//400)):
                    TF = n*math.log(n/(8*math.e*(m+n)))
                    nu = NU[m:m+n]
                    main = TF + lrmin - nu*math.log(d0); far = math.log(S1) - nu*math.log(D)
                    if np.any(far >= main): continue
                    low = main + np.log1p(-np.exp(far - main))
                    if np.all(low > U[m:m+n]): return (lam, Df, m, n)
    return None
def gmax(beta, mode, rig, b=B1):
    if not excluded(beta, 14.135, mode, rig, b): return None
    lo, hi = math.log(14.135), math.log(1e15)
    if excluded(beta, 1e15, mode, rig, b): return 1e15
    for _ in range(25):
        mid = (lo+hi)/2
        if excluded(beta, math.exp(mid), mode, rig, b): lo = mid
        else: hi = mid
    return math.exp(lo)
def yneeded(beta, g0, rig=False):
    lo, hi = math.log(10), math.log(10)*3000
    if not excluded(beta, g0, 'hyp', rig, hi): return None
    for _ in range(30):
        mid = (lo+hi)/2
        if excluded(beta, g0, 'hyp', rig, mid): hi = mid
        else: lo = mid
    return hi/math.log(10)
if __name__ == '__main__':
    print(f"S1 <= {S1:.5f}; C = {CRH:.7f}")
    what = sys.argv[1]
    if what == 'region':
        mode = sys.argv[2]; rig = sys.argv[3] == 'rig'
        if mode == 'window':
            a = math.log(1e6); C0 = float(env[int(a/W):].max()); mode = (a, C0)
            print(f"window [1e6,1e8], C0 = {C0:.6f}")
        for beta in [0.999, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55]:
            g = gmax(beta, mode, rig)
            print(f"  beta>={beta}: no zero with |gamma| <= {g if g is None else f'{g:.4g}'}", flush=True)
    else:
        rig = sys.argv[2] == 'rig'
        for g0 in [float(x) for x in sys.argv[3].split(',')]:
            row = []
            for beta in [0.51, 0.6, 0.75, 0.9, 0.99]:
                lys = math.log10(CRH*abs(complex(beta, g0)*(1-complex(beta, g0))))/(beta-0.5)
                yn = yneeded(beta, g0, rig)
                row.append(f"b={beta}: tube 10^{lys:.1f} / one-sided 10^{yn if yn is None else f'{yn:.1f}'}")
            print(f"t={g0:.4g}: " + "; ".join(row), flush=True)
