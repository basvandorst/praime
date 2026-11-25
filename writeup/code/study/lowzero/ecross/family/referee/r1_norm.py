"""Referee check 1: harmonic weights via PARI's Petersson norm (mfpetersson, normalised by 1/[SL2Z:Gamma0(N)]),
the Petersson sum rule sum_f omega_f lam_f(n) = Delta(1,n) with an independent Kloosterman code, the identity
omega_f = 2 pi^2/((k-1) N L(1,Sym^2 f)) (L from lfunsympow for rational k=2 forms), root numbers eps_f = -i^k lam_f(N) sqrt N
against lfunrootres for k = 2,4,6,8,10,14 (both signs)."""
import cypari2, math, time
import numpy as np
from scipy.special import jv
pari = cypari2.Pari(); pari.allocatemem(600*10**6); pari.set_real_precision(38)

def kl(m, n, c):
    if c == 1: return 1.0
    x = np.arange(1, c, dtype=np.int64); x = x[np.gcd(x, c) == 1]
    inv = np.array([pow(int(a), -1, c) for a in x], dtype=np.int64)
    return float(np.cos(2*np.pi*((m*x + n*inv) % c)/c).sum())

def Delta(k, N, m, n, C):
    s = 0.0
    for c in range(N, C+1, N):
        s += kl(m, n, c)/c*jv(k-1, 4*math.pi*math.sqrt(m*n)/c)
    return (1.0 if m == n else 0.0) + 2*math.pi*(-1)**(k//2)*s

def forms(N, k, P):
    pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
    out = []
    for j in range(1, int(pari('#B'))+1):
        co = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, P))
        pet = pari('mfpetersson(mfsymbol(mf,B[%d]))' % j)
        if str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) == 't_VEC':
            embs = [(list(co[i]), float(pet[i, i].real())) for i in range(len(co))]
        else:
            embs = [(list(co), float(pet.real()))]
        for (c, pv) in embs:
            a = np.array([float(complex(x).real) for x in c])
            lam = np.zeros(len(a)); lam[1:] = a[1:]/np.arange(1, len(a))**((k-1)/2)
            om = math.gamma(k-1)/((4*math.pi)**(k-1)*(N+1)*pv)
            out.append(dict(j=j, lam=lam, om=om))
    return out

for (N, k, C) in []:
    t0 = time.time()
    fs = forms(N, k, 12)
    print('N=%d k=%d: %d forms, omega_f (mfpetersson) = %s' % (N, k, len(fs), [round(f['om'], 6) for f in fs]))
    for n in [1, 2, 3, 5, 6]:
        lhs = sum(f['om']*f['lam'][n] for f in fs)
        rhs = Delta(k, N, 1, n, C)
        print('   n=%d  sum omega lam(n) = %.9f   Delta(1,n) Kloosterman (c<=%d) = %.9f  diff %.1e' % (n, lhs, C, rhs, lhs-rhs))
    print('   %.0fs' % (time.time()-t0))

# omega vs L(1,Sym^2) for rational weight-2 forms (elliptic curves), lfunsympow independent of the note's lfuncreate
for (lab, N, ai) in [('11a1', 11, [0,-1,1,-10,-20]), ('37a1', 37, [0,0,1,-1,0]), ('37b1', 37, [0,1,1,-23,-50]), ('43a1', 43, [0,1,1,0,0])]:
    E = pari.ellinit(ai); assert int(pari.ellglobalred(E)[0]) == N
    Ls = pari.lfunsympow(E, 2)
    L1 = float(pari.lfun(Ls, 2).real())   # motivic: s=2 <-> analytic s=1 for Sym^2 of weight 2
    fs = forms(N, 2, 40)
    ap = [int(pari.ellap(E, p)) for p in [2, 3, 5, 7]]
    for f in fs:
        if all(abs(f['lam'][p]*math.sqrt(p) - a) < 1e-6 for p, a in zip([2, 3, 5, 7], ap)):
            print('%s: omega(mfpetersson) = %.12f, 2pi^2/(N L(1,Sym^2 E)) with lfunsympow = %.12f' % (lab, f['om'], 2*math.pi**2/(N*L1)))

# root numbers: eps_f = -i^k lam_f(N) sqrt N, i.e. (-1)^{k/2+1} lam_f(N) sqrt N, against lfunrootres
for (N, k) in [(11, 6), (13, 6), (7, 8), (11, 10), (5, 14), (13, 10), (19, 6)]:
    pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
    res = []
    for j in range(1, int(pari('#B'))+1):
        co = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, N))
        Ls = pari('lfunmf(mf,B[%d])' % j)
        isvec = str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) == 't_VEC'
        embs = list(co) if isvec else [co]
        for i, c in enumerate(embs):
            lamN = float(complex(c[N]).real)/N**((k-1)/2)
            pred = (-1)**(k//2+1)*lamN*math.sqrt(N)
            L = Ls[i] if isvec else Ls
            r = complex(pari.lfunrootres(L)[2])
            res.append((round(pred, 6), round(r.real, 6)))
    print('N=%d k=%d (pred, PARI root number):' % (N, k), res)
