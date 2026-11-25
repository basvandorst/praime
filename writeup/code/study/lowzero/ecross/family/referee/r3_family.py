"""Referee check 3 (independent code): harmonic average of c_f with omega_f from PARI's Petersson norm (no L-function),
L(1,Sym^2 f) := 2pi^2/((k-1) N omega_f), B_f Euler product over p <= P; and a direct numerical test of the family main term
sum^h eps_f Z_f(H) = -(1/2pi) int H + errors, through the exact prime side of Proposition 1 (I_n by quadrature)."""
import cypari2, math, time, sys
import numpy as np
from scipy.special import loggamma, digamma
pari = cypari2.Pari(); pari.allocatemem(900*10**6); pari.set_real_precision(38)

def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

def forms(N, k, P):
    pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
    out = []
    for j in range(1, int(pari('#B'))+1):
        co = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, P))
        pet = pari('mfpetersson(mfsymbol(mf,B[%d]))' % j)
        isv = str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) == 't_VEC'
        embs = [(co[i], float(pet[i, i].real())) for i in range(len(co))] if isv else [(co, float(pet.real()))]
        for (c, pv) in embs:
            a = np.array([float(complex(x).real) for x in c])
            lam = np.zeros(len(a)); lam[1:] = a[1:]/np.arange(1, len(a))**((k-1)/2)
            om = math.gamma(k-1)/((4*math.pi)**(k-1)*(N+1)*pv)
            eps = int(round((-1)**(k//2+1)*lam[N]*math.sqrt(N)))
            out.append(dict(lam=lam, om=om, eps=eps))
    return out

def cfun(lam, N, k, om, P):
    pr = primes_upto(P); pr = pr[pr != N]
    l = lam[pr]; e = l*l - 1; p = pr.astype(float)
    num = 1 - e/p; den = 1 - e/p + e/p**2 - 1/p**3
    B = np.prod(np.sign(num))*math.exp(np.sum(np.log(np.abs(num)) - np.log(den)) - math.log1p(-1/N))
    LS = 2*math.pi**2/((k-1)*N*om)
    return B/LS, LS

which = sys.argv[1] if len(sys.argv) > 1 else 'c'
if which == 'c':
    for (N, k, P) in [(101, 4, 5000), (103, 4, 5000), (61, 4, 5000), (59, 4, 5000)]:
        t0 = time.time()
        fs = forms(N, k, P)
        om = np.array([f['om'] for f in fs]); eps = np.array([f['eps'] for f in fs])
        cs = np.array([cfun(f['lam'], N, k, f['om'], P)[0] for f in fs])
        print('N=%d k=%d: %d forms, sum omega = %.5f, sum^h eps = %+.5f, sum^h c_f = %.4f, natural mean c_f = %.4f '
              '(limit (6/pi^2)/(1-1/N) = %.4f), sum^h c_f^2 = %.3f, c_f range [%.3f, %.3f]  (%.0fs)'
              % (N, k, len(fs), om.sum(), np.sum(om*eps), np.sum(om*cs), cs.mean(), 6/math.pi**2/(1-1/N),
                 np.sum(om*cs*cs), cs.min(), cs.max(), time.time()-t0), flush=True)
else:
    N, k, P = 61, 4, 12000
    t0 = time.time()
    fs = forms(N, k, P)
    print('N=%d k=%d: %d forms, data %.0fs' % (N, k, len(fs), time.time()-t0), flush=True)
    om = np.array([f['om'] for f in fs]); eps = np.array([f['eps'] for f in fs])
    cs = np.array([cfun(f['lam'], N, k, f['om'], 5000)[0] for f in fs])
    lamM = np.array([f['lam'] for f in fs])           # forms x (P+1)
    pr = primes_upto(P)
    theta = lambda t: t/2*math.log(N) - t*math.log(2*math.pi) + np.imag(loggamma(k/2 + 1j*t))
    W = lambda u: np.where((u > 1) & (u < 2), np.exp(1 - 1/np.clip(1 - (2*u - 3)**2, 1e-300, None)), 0.0)
    for T in [4.0, 6.0, 8.0, 10.0]:
        tt = np.arange(T - 8, 2*T + 8, 0.004)
        ss = np.arange(T, 2*T, 0.004)
        Ws = W(ss/T)
        H = np.array([np.sum(Ws*np.exp(-(t - ss)**2))*0.004/math.sqrt(math.pi) for t in tt])
        intH = np.trapezoid(H, tt)
        ph = 2*theta(tt)
        dth = math.log(N) - 2*math.log(2*math.pi) + 2*np.real(digamma(k/2 + 1j*tt[(tt > T-3) & (tt < 2*T+3)]))
        nmin, nmax = math.exp(dth.min() - 2.5), math.exp(dth.max() + 2.5)
        assert nmax <= P, nmax
        tot = 0j; npp = 0
        for p in pr:
            p = int(p); q = p; j = 1
            while q <= nmax:
                if q >= nmin:
                    if p == N:
                        Lam = lamM[:, N]**j*math.log(N)
                    else:
                        lj = lamM[:, q]; lj2 = lamM[:, q//(p*p)] if j >= 2 else 0.0
                        Lam = (lj - lj2)*math.log(p)
                    cbar = np.sum(om*eps*Lam)
                    In = np.trapezoid(H*np.exp(1j*(ph - tt*math.log(q))), tt)
                    tot += cbar*q**-0.5*In; npp += 1
                q *= p; j += 1
        fam = -tot/(2*math.pi)
        print('T=%4.1f: n in [%.0f, %.0f] (%d prime powers); sum^h eps_f Z_f(H) = %.4f%+.4fi; -int H/2pi = %.4f; '
              'ratio %.4f; -(sum^h c_f) int H/2pi = %.4f  (sum^h c_f = %.4f)' % (T, nmin, nmax, npp, fam.real, fam.imag,
              -intH/(2*math.pi), fam.real/(-intH/(2*math.pi)), -np.sum(om*cs)*intH/(2*math.pi), np.sum(om*cs)), flush=True)
