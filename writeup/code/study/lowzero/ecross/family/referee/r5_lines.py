"""Referee check 5: Theorem D numerically at N = 61, k = 4: sum^h eps_f lam_f(l) Z_f(H; r) through the exact prime side of
Proposition 1 with the twist (I_n(r) by quadrature), against -(int H/2pi) A(l, r), A = mu(m) c_m(l)/(phi(m) sqrt l) if r = l/m^2."""
import cypari2, math, time
import numpy as np
from scipy.special import loggamma, digamma
pari = cypari2.Pari(); pari.allocatemem(900*10**6); pari.set_real_precision(38)
def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
N, k, P = 61, 4, 20000
t0 = time.time()
pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
lamM, om, eps = [], [], []
for j in range(1, int(pari('#B'))+1):
    co = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, P))
    pet = pari('mfpetersson(mfsymbol(mf,B[%d]))' % j)
    isv = str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) == 't_VEC'
    for i in (range(len(co)) if isv else [None]):
        c = co[i] if isv else co
        a = np.array([float(complex(x).real) for x in c])
        lam = np.zeros(len(a)); lam[1:] = a[1:]/np.arange(1, len(a))**((k-1)/2)
        pv = float(pet[i, i].real()) if isv else float(pet.real())
        lamM.append(lam); om.append(math.gamma(k-1)/((4*math.pi)**(k-1)*(N+1)*pv))
        eps.append(int(round((-1)**(k//2+1)*lam[N]*math.sqrt(N))))
lamM = np.array(lamM); om = np.array(om); eps = np.array(eps)
print('data: %d forms, %.0fs' % (len(om), time.time()-t0), flush=True)
pr = primes_upto(P)
theta = lambda t: t/2*math.log(N) - t*math.log(2*math.pi) + np.imag(loggamma(k/2 + 1j*t))
W = lambda u: np.where((u > 1) & (u < 2), np.exp(1 - 1/np.clip(1 - (2*u - 3)**2, 1e-300, None)), 0.0)
T = 6.0
tt = np.arange(T - 8, 2*T + 8, 0.004); ss = np.arange(T, 2*T, 0.004); Ws = W(ss/T)
H = np.array([np.sum(Ws*np.exp(-(t - ss)**2))*0.004/math.sqrt(math.pi) for t in tt])
intH = np.trapezoid(H, tt); ph = 2*theta(tt)
dth = math.log(N) - 2*math.log(2*math.pi) + 2*np.real(digamma(k/2 + 1j*tt[(tt > T-3) & (tt < 2*T+3)]))
lo, hi = math.exp(dth.min() - 2.5), math.exp(dth.max() + 2.5)
cases = [(1, 1/9, 0.5), (2, 2.0, 1/math.sqrt(2)), (2, 0.5, -1/math.sqrt(2)), (4, 1.0, -0.5),
         (3, 3.0, 1/math.sqrt(3)), (1, 0.5, 0.0), (1, 2.0, 0.0), (2, 1.0, 0.0), (1, 1.5, 0.0)]
for (l, r, A) in cases:
    nmin, nmax = lo/r, hi/r
    if nmax > P:
        print("l=%d r=%.4f: skipped (n up to %.0f > P)" % (l, r, nmax)); continue
    wl = om*eps*lamM[:, l]
    tot = 0j
    for p in pr:
        p = int(p); q = p; j = 1
        while q <= nmax:
            if q >= nmin:
                if p == N: Lam = lamM[:, N]**j*math.log(N)
                else: Lam = (lamM[:, q] - (lamM[:, q//(p*p)] if j >= 2 else 0.0))*math.log(p)
                In = np.trapezoid(H*np.exp(1j*(ph - tt*math.log(q*r))), tt)
                tot += np.sum(wl*Lam)*q**-0.5*In
            q *= p; j += 1
    fam = -tot/(2*math.pi)
    print('l=%d r=%.4f: family %.4f%+.4fi ; theorem -(int H/2pi) A = %.4f  (A = %.4f)' % (l, r, fam.real, fam.imag,
          -intH/(2*math.pi)*A, A), flush=True)
