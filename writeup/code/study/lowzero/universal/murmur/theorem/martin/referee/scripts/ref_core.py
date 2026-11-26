# Referee's own implementation (k = 2 only): Martin Prop 5.5 traces, limit densities for P = 3.
import math
from fractions import Fraction
import cypari2
pari = cypari2.Pari(); pari.allocatemem(400*10**6)

_H = {}
def H(n):
    """Hurwitz class number H(n) (n>0), 0 if n = 1,2 mod 4; H(0) not used."""
    if n % 4 in (1, 2): return Fraction(0)
    if n not in _H:
        _H[n] = Fraction(str(pari.qfbhclassno(n)))
    return _H[n]

def kron(a, b): return int(pari.kronecker(a, b))

def xi(Delta, v):
    D = -Delta
    if Delta % (v*v) != 0:
        return Fraction(kron(Delta, v) - 1)
    h = H(D)
    return v*H(D//(v*v))/h - 1

def primefactors(n): return [int(p) for p in pari.factor(n)[0]]
def mu(n):
    f = pari.factor(n); 
    if n == 1: return 1
    if any(int(e) > 1 for e in f[1]): return 0
    return (-1)**len(f[0])

def tQ(N, Q, l):
    """l^(1-k/2) tr(T_l W_Q | S_2^new(N)), N squarefree, l prime not dividing N (Martin Prop 5.5, k = 2)."""
    M = N // Q
    Mp = primefactors(M) if M > 1 else []
    tot = Fraction(0)
    s = 0
    while s*s*Q < 4*l:
        Delta = s*s*Q*Q - 4*l*Q
        x = Fraction(1)
        for v in Mp: x *= xi(Delta, v)
        term = x*H(-Delta)
        tot += term if s == 0 else 2*term
        s += 1
    val = -tot/2 + mu(M)*(l+1)
    if N == 1: val -= 1
    return val

# ---- local factors (closed forms as stated in [ML]/[AL]); the referee re-derives them numerically elsewhere
def e_v(v, s):
    if s % v == 0: return v*v/(v*v-1.0)
    return 1 - 1.0/((v-1)**2*(v+1))
def e_new(v, s):
    if s % v == 0: return -v/(v+1.0)
    return -v*v/(v*v-1.0)
PR = [int(p) for p in pari.primes(20000)]
def nuprime(s):
    r = 1.0
    for v in PR: r *= e_v(v, s)
    return r
ZETA2 = math.pi**2/6
_nu = {}
def nup(s):
    if s not in _nu: _nu[s] = nuprime(s)
    return _nu[s]

def D_dens(y, Q0, P=3):
    """D^[Q0,M](y), k = 2, P = 3 (note's (1.1)); kappa(P) = P/phi(P)."""
    M = P // Q0
    kap = P/ (P-1)
    E0 = ZETA2
    for v in primefactors(Q0) if Q0 > 1 else []: E0 *= (1 - v**-2)
    for v in primefactors(M) if M > 1 else []: E0 *= -(v-1)/v
    tot = -E0*math.sqrt(y/M)        # (-1)^(k/2) = -1 for k = 2
    ssum = 0.0
    s = 1
    while s*s < 4*y*M:
        Es = nup(s)
        for v in primefactors(P): Es /= e_v(v, s)
        for v in (primefactors(M) if M > 1 else []): Es *= e_new(v, s)
        ssum += Es*math.sqrt(4*y*M - s*s)
        s += 1
    tot -= ssum/M
    return 12*kap/math.pi*tot + 12*kap*mu(M)*y

def M_eps(x, beta, e3, eq, P=3, n=4000):
    """eps_q * sum_{Q0|3} eps(Q0) int_1^beta D(x/t) t dt / int t dt (midpoint rule)."""
    num = 0.0; den = 0.0
    h = (beta-1)/n
    for i in range(n):
        t = 1 + (i+0.5)*h
        num += (D_dens(x/t, 3) * e3 + D_dens(x/t, 1)) * t
        den += t
    return eq*num/den
