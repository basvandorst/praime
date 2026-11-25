"""Referee check 2 (independent code): Delta data point and the numbers of point 6.
- zeros of L(Delta) to 420 by lfunzeros; completeness by sign changes of Z(t) (lfunhardy, step 0.02) and N(T) - theta(T)/pi;
- theta(t) = -t ln 2pi + Im lnGamma(6+it) (mpmath), one-point sums sharp/Hann on (20,400];
- <Delta,Delta> by mfpetersson, omega_Delta by an own Kloosterman series, L(1,Sym^2 Delta) = 2pi^2/(11 omega);
- B_Delta and the direct product to P = 10^6; c_Delta;
- point 6: classical main terms sum_c mu(c)^2/phi(c) Hhat(2 ln c)/Hhat(0); Hardy-Littlewood series sum_q mu(q) lam(q^2)/q."""
import cypari2, math, time
import numpy as np, mpmath as mp
from scipy.special import jv
pari = cypari2.Pari(); pari.allocatemem(700*10**6); pari.set_real_precision(38)
t0 = time.time()
pari('mf=mfinit([1,12],1); D=mfDelta(); L=lfunmf(mf,D);')
g = np.array([float(x) for x in pari('lfunzeros(L,420)')])
print('zeros: %d in [0,420], first %.6f, min gap %.4f (%.0fs)' % (len(g), g[0], np.min(np.diff(g)), time.time()-t0), flush=True)
t0 = time.time()
Z = np.array([float(x) for x in pari('Li=lfuninit(L,[425]); vector(21000,i,lfunhardy(Li,i*0.02))')])
tt = 0.02*np.arange(1, 21001)
sc = np.sum(np.sign(Z[1:]) != np.sign(Z[:-1]))
print('sign changes of Z on (0.02,420]: %d; lfunzeros count %d (%.0fs)' % (sc, len(g), time.time()-t0), flush=True)
mp.mp.dps = 25
th = lambda t: float(-t*mp.log(2*mp.pi) + mp.im(mp.loggamma(6 + 1j*t)))
ths = np.array([th(x) for x in g])
dev = [int(np.sum(g <= T)) - th(T)/math.pi for T in np.arange(30, 420, 1.0)]
print('N(T) - theta(T)/pi on T = 30..419: mean %.3f, min %.3f, max %.3f' % (np.mean(dev), np.min(dev), np.max(dev)))
for (T1, T2) in [(20, 400), (20, 140), (140, 260), (260, 400)]:
    m = (g > T1) & (g <= T2)
    sh = np.sum(np.exp(2j*ths[m]))
    hn = np.sum(np.sin(math.pi*(g[m]-T1)/(T2-T1))**2*np.exp(2j*ths[m]))
    print('[%d,%d]: %d zeros, sharp %.3f%+.3fi, Hann %.3f%+.3fi' % (T1, T2, m.sum(), sh.real, sh.imag, hn.real, hn.imag))
# Petersson norm and omega
pet = float(pari('mfpetersson(mfsymbol(mf,D))'))
om_pet = math.gamma(11)/((4*math.pi)**11*pet)
def kl(m, n, c):
    if c == 1: return 1.0
    x = np.arange(1, c, dtype=np.int64); x = x[np.gcd(x, c) == 1]
    inv = np.array([pow(int(a), -1, c) for a in x], dtype=np.int64)
    return float(np.cos(2*np.pi*((m*x + n*inv) % c)/c).sum())
om_kl = 1 + 2*math.pi*sum(kl(1, 1, c)/c*jv(11, 4*math.pi/c) for c in range(1, 60))
LS = 2*math.pi**2/(11*om_pet)
print('<Delta,Delta> = %.15e; omega via norm %.13f, via Kloosterman series %.13f; L(1,Sym^2) = %.12f; 1/omega = %.6f'
      % (pet, om_pet, om_kl, LS, 1/om_pet))
# Euler products, P = 10^6 (in GP)
t0 = time.time()
res = pari('v=mfcoefs(D,10^6); lb=0.; ld=0.; sg=1; forprime(p=2,10^6, l=v[p+1]/p^(11/2); e=l^2-1; '
           'lb+=log(1-e/p)-log(1-e/p+e/p^2-1/p^3); ld+=log(abs(1-e/p)); if(1-e/p<0,sg=-sg)); [exp(lb), sg*exp(ld)]')
B, direct = float(res[0]), float(res[1])
print('B_Delta (p<=1e6) = %.10f; c_Delta = B/L = %.6f; direct prod_{p<=1e6}(1-lam(p^2)/p) = %.6f; omega*c = %.4f (%.0fs)'
      % (B, B/LS, direct, om_pet*B/LS, time.time()-t0), flush=True)
cD = B/LS
for (T1, T2) in [(20, 400)]:
    print('law sharp -c T/2pi = %.3f, Hann %.3f; c=1 amplitude only: %.3f' % (-cD*(T2-T1)/(2*math.pi), -cD*(T2-T1)/(4*math.pi),
                                                                        -(T2-T1)/(2*math.pi*om_pet)))
# lambda(p) for p <= 10^6 into numpy (primes only) for the HL series
lamp = np.array([float(x) for x in pari('vector(#primes([2,10^6]),i,my(p=prime(i)); v[p+1]/p^(11/2))')])
P = np.array([int(x) for x in pari('primes([2,10^6])')])
Q = 10**6
gq = np.ones(Q+1)
for p, l in zip(P, lamp):
    gq[p::p] *= -(l*l - 1)
    if p*p <= Q: gq[p*p::p*p] = 0.0
cum = np.cumsum(gq[1:]/np.arange(1, Q+1))
ep = np.cumprod(1 - (lamp**2 - 1)/P)
for q in [1, 2, 3, 5, 10, 30, 100, 1000, 10**4, 10**5, 10**6]:
    print('  HL partial sum sum_{q<=%7d} mu(q) lam(q^2)/q = %.4f   Euler partial prod_{p<=%d}(1-lam(p^2)/p) = %.4f'
          % (q, cum[q-1], q, ep[np.searchsorted(P, q, side='right')-1] if q >= 2 else 1.0))
# point 6 (a): classical PNT-in-AP main terms of the Petersson moduli: sum_c mu(c)^2/phi(c) Hhat(2 ln c)/Hhat(0), Hann [20,400]
T1, T2 = 20.0, 400.0
s = np.linspace(T1, T2, 400001); w = np.sin(np.pi*(s-T1)/(T2-T1))**2
def Hhat(xi): return np.trapezoid(w*np.exp(-1j*xi*s), s)
def phi(n):
    r, m, p = n, n, 2
    while p*p <= m:
        if m % p == 0:
            while m % p == 0: m //= p
            r -= r//p
        p += 1
    if m > 1: r -= r//m
    return r
def sqfree(n):
    p = 2
    while p*p <= n:
        if n % (p*p) == 0: return False
        p += 1
    return True
H0 = Hhat(0.0).real
tot = 0j; big = 0.0
for c in range(1, 2001):
    if sqfree(c):
        v = Hhat(2*math.log(c))/H0/phi(c)
        tot += v
        if c >= 2: big = max(big, abs(v))
print('point 6(a): sum_{c<=2000} mu(c)^2/phi(c) Hhat(2 ln c)/Hhat(0) = %.8f%+.8fi (max |c>=2 term| %.2e): classical amplitude '
      '= %.6f x (1/omega) = %.6f, against c_Delta = %.6f' % (tot.real, tot.imag, big, tot.real, tot.real/om_pet, cD))
