# check: int_1^inf H(y) y^(-s-1) dy = -zeta'/zeta(s)/(s(s-1)) - 1/(s-1)^2 + (1+g)/(s-1) - ln(2pi)/s + sum_k 1/(2k(2k+1)(s+2k))
import mpmath as mp, numpy as np, math
mp.mp.dps = 25
N = 10**7
s_ = np.ones(N+1, bool); s_[:2] = False
for p in range(2, int(N**.5)+1):
    if s_[p]: s_[p*p::p] = False
P = np.nonzero(s_)[0]; del s_
pp = []
for p in P:
    q = int(p); lp = math.log(q); x = q
    while x <= N: pp.append((x, lp)); x *= q
    if q*q > N: pp.pop(); break
pp += [(int(p), math.log(int(p))) for p in P[P > int(N**.5)]]
pp.sort()
n = np.array([a for a,_ in pp], float); w = np.array([b for _,b in pp])
g = float(mp.euler)
def F_num(s):
    # piecewise exact integral of G y^(-s-1) on [1,N], plus -ln2pi*int + T part done analytically
    A = np.concatenate([[0.0], np.cumsum(w/n)]); Psi = np.concatenate([[0.0], np.cumsum(w)])
    a = np.concatenate([[1.0], n]); b = np.concatenate([n, [float(N)]])
    # int y^(-s) dy, int y^(-s-1), int y^(-s) ln y
    I1 = lambda u: u**(1-s)/(1-s); I2 = lambda u: -u**(-s)/s
    I3 = lambda u: u**(1-s)*(np.log(u)/(1-s) - 1/(1-s)**2)
    tot = np.sum((A+1+g)*(I1(b)-I1(a)) - Psi*(I2(b)-I2(a)) - (I3(b)-I3(a)))
    T = float(mp.nsum(lambda k: 1/(2*k*(2*k+1)*(s+2*k)), [1, mp.inf]))
    return tot - math.log(2*math.pi)*(1/s - N**(-s)/s) + T
def F_th(s):
    s = mp.mpf(s)
    zz = -mp.diff(mp.zeta, s)/mp.zeta(s)
    return zz/(s*(s-1)) - 1/(s-1)**2 + (1+mp.euler)/(s-1) - mp.log(2*mp.pi)/s + mp.nsum(lambda k: 1/(2*k*(2*k+1)*(s+2*k)), [1, mp.inf])
for s in (2.0, 3.0, 1.5):
    print(s, F_num(s), float(F_th(s)), "tail bound", 0.0462*N**(0.5-s)/(s-0.5))
