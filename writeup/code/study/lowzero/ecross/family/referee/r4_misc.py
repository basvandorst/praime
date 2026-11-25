"""Referee check 4: (i) S(Na,b;N^j) = 0 for j >= 2, (b,N) = 1; (ii) |c_chi| = |tau(chi)|^2/phi(m) for K(a) = S(Na,1;Nm);
(iii) Bessel-Mellin (*) by Abel damping e^{-eps y} (independent of quadosc); (iv) Delta: size of the Petersson c-terms against
their trivial and random (root-sum-square) scales; (v) Hardy-Littlewood: sum_n c_q(n) lam(n) w_n against (phi(q)/q) lam(q^2) sum_n lam(n) w_n."""
import math, cmath, numpy as np, mpmath as mp, cypari2
from scipy.special import jv, loggamma, digamma
TP = 2*math.pi
def units(c):
    x = np.arange(1, c, dtype=np.int64); x = x[np.gcd(x, c) == 1]
    return x, np.array([pow(int(a), -1, c) for a in x], dtype=np.int64)
def S(m, n, c):
    if c == 1: return 1.0
    x, inv = units(c); return float(np.cos(TP*((m*x + n*inv) % c)/c).sum())
# (i)
w = 0.0
for N in [3, 5, 7, 11]:
    for j in [2, 3]:
        for a in range(1, 12):
            for b in [1, 2, 4]:
                if b % N: w = max(w, abs(S(N*a, b, N**j)))
print('(i) max |S(Na,b;N^j)|, j=2,3, (b,N)=1: %.1e' % w)
# (ii) characters mod m via discrete logs for prime m; generic m via orthogonality on the group (brute force)
def chars(m):
    import itertools
    pari = cypari2.Pari()
    G = pari.znstar(m)
    cyc = [int(c) for c in G[1]]; gens = [int(pari.lift(g)) for g in G[2]]
    elog = {}
    for e in itertools.product(*[range(c) for c in cyc]):
        a = 1
        for g, ei in zip(gens, e): a = a*pow(g, ei, m) % m
        elog[a] = e
    x, _ = units(m)
    out = []
    for fexp in itertools.product(*[range(c) for c in cyc]):
        out.append({int(a): cmath.exp(2j*math.pi*sum(fexp[i]*elog[int(a)][i]/cyc[i] for i in range(len(cyc)))) for a in x})
    return x, out
worst = 0.0
for (N, m) in [(7, 9), (7, 11), (13, 12), (11, 15), (5, 16)]:
    x, chs = chars(m)
    K = {int(a): S(N*int(a), 1, N*m) for a in x}
    for ch in chs:
        chat = sum(K[a]*ch[a].conjugate() for a in K)/len(x)
        tau = sum(ch[a]*cmath.exp(2j*math.pi*a/m) for a in K)
        worst = max(worst, abs(abs(chat) - abs(tau)**2/len(x)))
print('(ii) max ||c_chi| - |tau(chi)|^2/phi(m)| over 5 moduli, all chi: %.1e' % worst)
# (iii)
mp.mp.dps = 15
for (k, N, m, t) in [(2, 37, 1, 9.0), (12, 1, 1, 15.0), (12, 1, 3, 7.0)]:
    a = 4*mp.pi/(mp.sqrt(N)*m); s = mp.mpc(1, -2*t); nu = k-1
    vals = []
    for eps in [0.04, 0.02, 0.01]:
        f = lambda y: y**(s-1)*mp.besselj(nu, y)*mp.exp(-eps*y)
        vals.append(2*a**(-s)*mp.quad(f, [0, 0.5, 2] + list(mp.linspace(4, 60/eps, int(60/eps/8)+1))))
    rich = (8*vals[2] - 6*vals[1] + vals[0])/3      # Richardson for O(eps) + O(eps^2)
    th = t/2*mp.log(N*m*m) - t*mp.log(2*mp.pi) + mp.im(mp.loggamma(k/2 + 1j*t))
    ex = mp.sqrt(N)*m/(2*mp.pi)*mp.exp(-2j*th)
    print('(iii) k=%d N=%d m=%d t=%g: Abel-damped (Richardson) %s, (sqrt N m/2pi) e^{-2i theta} %s, rel diff %.1e'
          % (k, N, m, t, mp.nstr(rich, 10), mp.nstr(ex, 10), float(abs(rich-ex)/abs(ex))))
# (iv), (v): Delta, Hann [20,400], stationary weights at the exact point
pari = cypari2.Pari(); pari.allocatemem(300*10**6)
X = 5400
tau = np.array([float(int(v)) for v in pari('mfcoefs(mfDelta(),%d)' % X)])
n = np.arange(X+1, dtype=float); lam = np.zeros(X+1); lam[1:] = tau[1:]/n[1:]**5.5
nn = np.arange(2, X+1, dtype=float)
t = TP*np.sqrt(nn)
f = lambda t: 2*np.real(digamma(6 + 1j*t)) - 2*math.log(TP) - np.log(nn)
for _ in range(40):
    h = 1e-5; d = (f(t+h) - f(t-h))/(2*h); t = t - f(t)/d
phi2 = (f(t+1e-5) - f(t-1e-5))/2e-5
phase = 2*(-t*math.log(TP) + np.imag(loggamma(6 + 1j*t))) - t*np.log(nn)
H = np.where((t > 20) & (t < 400), np.sin(math.pi*(t-20)/380)**2, 0.0)
wn = np.zeros(X+1, complex)
wn[2:] = -(1/TP)*nn**-0.5*np.sqrt(TP/phi2)*np.exp(1j*(phase + math.pi/4))*H
R1 = np.sum(lam*wn)
print('(v) R_1 = sum_n lam(n) w_n = %.3f%+.3fi (Voronoi m=1 main term -(1/2pi) int H = %.3f)' % (R1.real, R1.imag, -190/TP))
def cq(q, nv):
    return np.array([sum(d*int(np.sign(1))*(1 if True else 0)*mob(q//d) for d in range(1, q+1) if q % d == 0 and v % d == 0) for v in nv])
def mob(n):
    r, m, p = 1, n, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
def phi(n): return sum(1 for a in range(1, n+1) if math.gcd(a, n) == 1)
idx = np.arange(X+1)
for q in [2, 3, 5, 6, 7, 10]:
    c = np.zeros(X+1)
    for d in range(1, q+1):
        if q % d == 0: c[idx % d == 0] += d*mob(q//d)
    Rq = np.sum(c*lam*wn)
    pred = phi(q)/q*lam[q*q]*R1
    print('   q=%2d: sum_n c_q(n) lam(n) w_n = %8.3f%+.3fi ; (phi(q)/q) lam(q^2) R_1 = %8.3f%+.3fi' % (q, Rq.real, Rq.imag, pred.real, pred.imag))
# (iv) Petersson c-terms over primes
om = 2.8402873751675
pr = np.array([p for p in range(2, X+1) if all(p % r for r in range(2, int(p**0.5)+1))])
for c in [2, 3, 5, 10, 20, 30, 40, 60]:
    if c == 1: Sc = np.ones(len(pr))
    else:
        Sa = np.array([S(a, 1, c) for a in range(c)]); Sc = Sa[pr % c]
    term = TP/om*Sc/c*jv(11, 4*math.pi*np.sqrt(pr)/c)*np.log(pr)*wn[pr]
    print('(iv) c=%2d: c-term %7.3f%+.3fi ; trivial sum|.| %7.2f ; root-sum-square %6.3f ; |c-term|/rss %.2f'
          % (c, term.sum().real, term.sum().imag, np.abs(term).sum(), math.sqrt(np.sum(np.abs(term)**2)),
             abs(term.sum())/math.sqrt(np.sum(np.abs(term)**2))))
