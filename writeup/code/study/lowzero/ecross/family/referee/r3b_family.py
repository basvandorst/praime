"""Referee check 3b: harmonic weights omega_f from the Petersson formula itself (least squares on
sum_f omega_f lam_f(n) = Delta_{k,N}(1,n), n <= 60, own Kloosterman code, k = 4), no L-function and no mfpetersson;
then L(1,Sym^2 f) = 2pi^2/((k-1) N omega_f), c_f = B_f/L(1,Sym^2 f), sum^h c_f. Cross-check with mfpetersson where it fits."""
import cypari2, math, time, sys
import numpy as np
from scipy.special import jv
pari = cypari2.Pari(); pari.allocatemem(800*10**6); pari.set_real_precision(38)

def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

def coeffs(N, k, P, petersson=False):
    pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
    out = []
    for j in range(1, int(pari('#B'))+1):
        co = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, P))
        isv = str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) == 't_VEC'
        if petersson:
            pet = pari('mfpetersson(mfsymbol(mf,B[%d]))' % j)
        embs = list(range(len(co))) if isv else [None]
        for i in embs:
            c = co[i] if isv else co
            a = np.array([float(complex(x).real) for x in c])
            lam = np.zeros(len(a)); lam[1:] = a[1:]/np.arange(1, len(a))**((k-1)/2)
            omp = None
            if petersson:
                pv = float(pet[i, i].real()) if isv else float(pet.real())
                omp = math.gamma(k-1)/((4*math.pi)**(k-1)*(N+1)*pv)
            out.append(dict(lam=lam, omp=omp, eps=int(round((-1)**(k//2+1)*lam[N]*math.sqrt(N)))))
    return out

def Delta_vec(k, N, nmax, M):
    ns = np.arange(1, nmax+1)
    s = np.zeros(nmax)
    for c in range(N, N*M+1, N):
        x = np.arange(1, c, dtype=np.int64); x = x[np.gcd(x, c) == 1]
        inv = np.array([pow(int(a), -1, c) for a in x], dtype=np.int64)
        S = np.cos(2*np.pi*((x[None, :] + (ns[:, None] % c)*inv[None, :]) % c)/c).sum(axis=1)
        s += S/c*jv(k-1, 4*math.pi*np.sqrt(ns)/c)
    return (ns == 1).astype(float) + 2*math.pi*(-1)**(k//2)*s

def cB(lam, N, P):
    pr = primes_upto(P); pr = pr[pr != N]
    l = lam[pr]; e = l*l - 1; p = pr.astype(float)
    num = 1 - e/p; den = 1 - e/p + e/p**2 - 1/p**3
    return np.prod(np.sign(num))*math.exp(np.sum(np.log(np.abs(num)) - np.log(den)) - math.log1p(-1/N))

k, P, nmax = 4, 5000, 60
for N in [int(a) for a in sys.argv[1:]]:
    t0 = time.time()
    fs = coeffs(N, k, P, petersson=(N <= 61))
    A = np.array([[f['lam'][n] for f in fs] for n in range(1, nmax+1)])
    rhs = Delta_vec(k, N, nmax, 300 if N < 80 else 200)
    om, res, rank, sv = np.linalg.lstsq(A, rhs, rcond=None)
    resid = np.max(np.abs(A @ om - rhs))
    eps = np.array([f['eps'] for f in fs])
    Bs = np.array([cB(f['lam'], N, P) for f in fs])
    LS = 2*math.pi**2/((k-1)*N*om)
    cs = Bs/LS
    line = ('N=%d k=%d: %d forms, rank %d, cond %.1e, max residual %.1e; sum omega = %.5f; sum^h eps = %+.5f; '
            'sum^h c_f = %.4f; natural mean c_f = %.4f (limit %.4f); sum^h c_f^2 = %.3f; c_f in [%.3f, %.3f]'
            % (N, k, len(fs), rank, sv[0]/sv[-1], resid, om.sum(), np.sum(om*eps), np.sum(om*cs), cs.mean(),
               6/math.pi**2/(1-1/N), np.sum(om*cs*cs), cs.min(), cs.max()))
    if fs[0]['omp'] is not None:
        omp = np.array([f['omp'] for f in fs])
        line += '; max |omega(lstsq) - omega(mfpetersson)| = %.1e' % np.max(np.abs(om - omp))
    print(line + '  (%.0fs)' % (time.time()-t0), flush=True)
