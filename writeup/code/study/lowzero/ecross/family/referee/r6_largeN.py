"""Referee check 6: Theorems B and D at a large prime level (N = 10007, k = 4), where the asymptotic regime is reachable.
Family coefficients from the root-number Petersson formula (Lemma 2, a theorem: Lemma 1 + Hecke + root number, each checked
in r1): G(n) := sum^h eps_f lam_f(n) = 2pi N^{-1/2} sum_{m<=M} S(n, Nbar; m)/m J_{k-1}(4pi sqrt(n/N)/m)
(S(Nn,1;Nm) = -S(n,Nbar;m), (m,N) = 1; the m = 0 mod N terms vanish identically). Then the exact prime side of Proposition 1:
sum^h eps_f lam_f(l) Z_f(H; r) = -(1/2pi) sum_n [sum^h eps_f lam_f(l) Lambda_f(n)] n^{-1/2} I_n(r) (n = p, p^2), I_n(r) by quadrature."""
import math, time, sys
import numpy as np
from scipy.special import jv, loggamma, digamma
N, k, T, M = 10007, 4, float(sys.argv[1]) if len(sys.argv) > 1 else 6.0, 1000
def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
t0 = time.time()
tables = {}
for m in range(1, M+1):
    if m == 1: tables[1] = np.array([1.0]); continue
    x = np.arange(1, m, dtype=np.int64); x = x[np.gcd(x, m) == 1]
    inv = np.array([pow(int(a), -1, m) for a in x], dtype=np.int64)
    Nb = pow(N, -1, m)
    a = np.arange(m, dtype=np.int64)
    tables[m] = np.cos(2*np.pi*((np.outer(a, x) + Nb*inv[None, :]) % m)/m).sum(axis=1)
print('Kloosterman tables m <= %d: %.0fs' % (M, time.time()-t0), flush=True)
def G(n):
    n = np.asarray(n, dtype=np.int64); s = np.zeros(len(n))
    y0 = 4*math.pi*np.sqrt(n/N)
    for m in range(1, M+1):
        s += tables[m][n % m]/m*jv(k-1, y0/m)
    return 2*math.pi/math.sqrt(N)*s
theta = lambda t: t/2*math.log(N) - t*math.log(2*math.pi) + np.imag(loggamma(k/2 + 1j*t))
W = lambda u: np.where((u > 1) & (u < 2), np.exp(1 - 1/np.clip(1 - (2*u - 3)**2, 1e-300, None)), 0.0)
tt = np.arange(T - 8, 2*T + 8, 0.01); ss = np.arange(T, 2*T, 0.002); Ws = W(ss/T)
H = np.array([np.sum(Ws*np.exp(-(t - ss)**2))*0.002/math.sqrt(math.pi) for t in tt])
intH = np.trapezoid(H, tt); A0 = H*np.exp(2j*theta(tt))
dth = math.log(N) - 2*math.log(2*math.pi) + 2*np.real(digamma(k/2 + 1j*tt[(tt > T-3) & (tt < 2*T+3)]))
lo, hi = math.exp(dth.min() - 1.8), math.exp(dth.max() + 1.8)
print('T=%g: int H = %.4f, -int H/2pi = %.4f; n*r in [%.0f, %.0f]' % (T, intH, -intH/(2*math.pi), lo, hi), flush=True)
def Iq(nr):
    out = np.zeros(len(nr), complex)
    for i in range(0, len(nr), 400):
        out[i:i+400] = np.trapezoid(A0[None, :]*np.exp(-1j*np.outer(np.log(nr[i:i+400]), tt)), tt, axis=1)
    return out
def mu(n):
    r, m, p = 1, n, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
def ram(m, l): return sum(d*mu(m//d) for d in range(1, m+1) if m % d == 0 and l % d == 0)
def phi(m): return sum(1 for a in range(1, m+1) if math.gcd(a, m) == 1)
def Aval(l, r):
    m2 = l/r; m = round(math.sqrt(m2))
    if m >= 1 and abs(m*m - m2) < 1e-9: return mu(m)*ram(m, l)/(phi(m)*math.sqrt(l))
    return 0.0
cases = [(1, 1.0), (1, 0.25), (2, 2.0), (2, 0.5), (4, 1.0), (3, 3.0), (1, 0.5), (2, 1.0), (1, 2.0), (1, 1.5)]
if len(sys.argv) > 2: cases = [cases[int(i)] for i in sys.argv[2].split(',')]
for (l, r) in cases:
    t1 = time.time()
    pr = primes_upto(int(hi/r) + 1); pr = pr[(pr >= lo/r) & (np.gcd(pr, l) == 1)]
    cp = np.log(pr)*G(l*pr)                              # sum^h eps lam(l) Lambda_f(p), p not dividing l
    tot = np.sum(cp*pr**-0.5*Iq(pr*r))
    p2 = primes_upto(int(math.sqrt(hi/r)) + 1); p2 = p2[(p2*p2 >= lo/r) & (np.gcd(p2, l) == 1)]
    if len(p2):
        c2 = np.log(p2)*(G(l*p2*p2) - G(np.array([l]))[0])  # Lambda_f(p^2) = (lam(p^2) - 1) ln p
        tot += np.sum(c2/p2*Iq(p2.astype(float)**2*r))
    fam = -tot/(2*math.pi); A = Aval(l, r)
    Hh = lambda xi: np.trapezoid(H*np.exp(-1j*xi*tt), tt)
    lines = -sum(mu(m)*ram(m, l)/(phi(m)*math.sqrt(l))*Hh(math.log(r*m*m/l)) for m in range(1, 61) if mu(m) != 0)/(2*math.pi)
    print('l=%d r=%.4f: %6d primes; family %.4f%+.4fi ; theorem -(int H/2pi) A = %.4f (A = %+.4f); finite-T line sum '
          '%.4f%+.4fi; family - lines = %.4f%+.4fi  [%.0fs]' % (l, r, len(pr), fam.real, fam.imag, -intH/(2*math.pi)*A, A,
          lines.real, lines.imag, (fam-lines).real, (fam-lines).imag, time.time()-t1), flush=True)
