# Independent check of Part II, Theorem 6.1 / Table 1: a_n = n S_n + delta_{n1}/4,
# S_n = sum_c (1/c)[S(1,n;c) J2dot(4 pi sqrt n / c) - S(-1,n;c) I2dot(4 pi sqrt n / c)].
import numpy as np, mpmath as mp
from math import gcd
mp.mp.dps = 25
def Jd(w): return mp.diff(lambda nu: mp.besselj(nu, w), 2)
def Id(w): return mp.diff(lambda nu: mp.besseli(nu, w), 2)
# check Lemma B.5(a): -I2dot = K2 + 2/w I1 - 2/w^2 I0, and Lemma B.1 series
for w in [0.3, 1.7, 5.0, 4*mp.pi]:
    lhs = -Id(w); rhs = mp.besselk(2,w) + 2/w*mp.besseli(1,w) - 2/w**2*mp.besseli(0,w)
    ser = mp.nsum(lambda k: (w/2)**(2*k+2)/(mp.factorial(k)*mp.factorial(k+2))*(mp.log(w/2)-mp.digamma(k+3)), [0, mp.inf])
    print("w=%s  -I2dot=%s  B.5(a) rhs=%s  B.1 series I2dot=%s" % (mp.nstr(w,6), mp.nstr(lhs,15), mp.nstr(rhs,15), mp.nstr(ser,15)))
def kloost(m, n, c):
    if c == 1: return 1.0
    d = np.array([x for x in range(c) if gcd(x, c) == 1])
    dinv = np.array([pow(int(x), -1, c) for x in d])
    return float(np.sum(np.cos(2*np.pi*((m*d + n*dinv) % c)/c)))
CMAX = 1500
table = {1: 4581.370553, 2: 1025727.23492, 3: 62590223.68999, 4: 1968035704.4549, 5: 40671323483.6182}
for n in range(1, 6):
    z = 4*mp.pi*mp.sqrt(n)
    S = mp.mpf(0); T1 = None; Rabs = mp.mpf(0)
    for c in range(1, CMAX+1):
        w = z/c
        term = (kloost(1, n, c)*Jd(w) - kloost(-1, n, c)*Id(w))/c
        if c == 1: T1 = term
        else: Rabs += abs(term)
        S += term
    a = n*S + (mp.mpf(1)/4 if n == 1 else 0)
    abar = mp.mpf(n)**0.25*mp.e**z/(4*mp.sqrt(2)*mp.pi**2)
    print("n=%d a_n(c<=%d)=%s table=%s rel.diff=%.2e | a_n/abar=%s  (1-2/z)=%s  sum_{c>=2}|term|/T1=%.3e" % (
        n, CMAX, mp.nstr(a, 14), table[n], float((a-table[n])/table[n]), mp.nstr(a/abar, 8), mp.nstr(1-2/z, 8), float(Rabs/T1)))
# lambda_-(4 pi) and its pieces
z = 4*mp.pi; e = lambda w: mp.e**w/mp.sqrt(2*mp.pi*w)
k1 = 1 + mp.e/3; k2 = 2*mp.e/9; k3 = k1*mp.digamma(3) + k2
lam0 = (1-1/z)*(1-3/(8*z)-15/(32*z**2)) - mp.pi**3/(16*z**2)
pieces = [(mp.pi+2)*z/(8*e(z)), mp.sqrt(2)*mp.e**(-z/2), (mp.pi+2)*z**2/(8*e(z)), mp.sqrt(3)*z**2/4*mp.e**(-2*z/3), (k3 + z/2*(k1+k3))*z/(2*e(z))]
lm = lam0 - (mp.pi+1)*z/(4*e(z)) - sum(pieces)
print("lambda0(4pi)=%s  (pi+1)z/(4e)=%s pieces=%s  lambda_-(4pi)=%s" % (mp.nstr(lam0,8), mp.nstr((mp.pi+1)*z/(4*e(z)),6), [mp.nstr(p,6) for p in pieces], mp.nstr(lm,8)))
# lambda_- increasing: sample
import numpy as np
def lamminus(z):
    lam0 = (1-1/z)*(1-3/(8*z)-15/(32*z**2)) - mp.pi**3/(16*z**2)
    pieces = [(mp.pi+2)*z/(8*e(z)), mp.sqrt(2)*mp.e**(-z/2), (mp.pi+2)*z**2/(8*e(z)), mp.sqrt(3)*z**2/4*mp.e**(-2*z/3), (k3 + z/2*(k1+k3))*z/(2*e(z))]
    return lam0 - (mp.pi+1)*z/(4*e(z)) - sum(pieces)
zs = [4*mp.pi*(1+0.05*k) for k in range(0, 200)]
vals = [lamminus(x) for x in zs]
print("lambda_- monotone on sample [4pi, 4pi*11]:", all(vals[i] < vals[i+1] for i in range(len(vals)-1)), " min:", mp.nstr(min(vals),8))
# Bessel lemma spot checks on grids
bad = []
for w in np.linspace(0.05, 2.0, 40):
    w = mp.mpf(w); v = w**2/4; L = mp.log(2/w)
    B = w**2/8*((L+mp.digamma(3))*(1+v/3*mp.e**v) + v/9*(1+v)*mp.e**v)
    Bc = w**2/8*(k1*(L+mp.digamma(3)) + k2)
    if not (abs(Id(w)) <= B <= Bc and abs(Jd(w)) <= B and -Id(w) >= w**2/8*(L+mp.digamma(3))): bad.append(('B2', float(w)))
for w in np.linspace(0.1, 40, 120):
    w = mp.mpf(w)
    if abs(Jd(w)) > (mp.pi+1)/2: bad.append(('B3a', float(w)))
    if w >= 1:
        lo = 2/w*e(w)*((1-1/w)*(1-3/(8*w)-15/(32*w**2)) - mp.pi**3/(16*w**2)); hi = 2/w**2 + 2/w*(1-1/w)*e(w)
        if not (lo <= -Id(w) <= hi): bad.append(('B5c', float(w)))
    if w >= 2 and abs(Id(w)) > mp.mpf(1)/2 + 2/w*e(w): bad.append(('B5d', float(w)))
print("Bessel lemma violations on grids:", bad)
