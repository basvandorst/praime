# review 162: own evaluation of the registered limit density (definitions from derivation.txt, own code), on several binnings,
# with the degree-2 kernel comparator and the conditional models (observed local types at l <= L0, L0 = 5, 7, 13, 29, 47).
# M1(y) = int_0^y m = sum_poles r Kt(w0) y^(1-w0)/(1-w0) + (1/pi) Re int_0^T D(c+it) Kt(c+it) y^(1-w)/(1-w) dt - sqrt(y)   (c = 0.3)
import sys, math, time, json, numpy as np
from scipy.special import loggamma, zeta as szeta
R = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/'
T = float(sys.argv[1]) if len(sys.argv) > 1 else 1000.0; DT = 0.01; C0 = 0.3
t0 = time.time()

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for q in range(2, int(n**0.5) + 1):
        if s[q]: s[q * q::q] = False
    return np.nonzero(s)[0]

# --- Riemann zeta for complex arrays, Euler-Maclaurin (own implementation), N >= |s|/pi + 40, 16 Bernoulli terms
B2K = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330, 854513/138, -236364091/2730, 8553103/6,
       -23749461029/870, 8615841276005/14322, -7709321041217/510]
def zeta_c(s):
    s = np.asarray(s, complex); out = np.empty_like(s); order = np.argsort(np.abs(s.imag))
    for a in range(0, len(s), 2000):
        ii = order[a:a + 2000]; ss = s[ii]; N = int(np.abs(ss).max() / math.pi) + 40
        acc = np.zeros(len(ss), complex); n = np.arange(1, N, dtype=float); ln = np.log(n)
        for b in range(0, len(n), 400):
            acc += np.exp(-np.outer(ss, ln[b:b + 400])).sum(1)
        NS = np.exp(-ss * math.log(N)); acc += NS / 2 + N * NS / (ss - 1)
        poch = ss.copy(); pw = NS / N
        for k, B in enumerate(B2K):
            acc += B / math.factorial(2 * k + 2) * poch * pw
            poch = poch * (ss + 2 * k + 1) * (ss + 2 * k + 2); pw = pw / N**2
        out[ii] = acc
    return out
# self-check of zeta_c against scipy (real) and a known value
assert abs(zeta_c(np.array([2.0 + 0j]))[0] - math.pi**2 / 6) < 1e-12
z14 = zeta_c(np.array([0.5 + 14.134725141734693j]))[0]; assert abs(z14) < 1e-8, z14

LNPI = math.log(math.pi)
def Kt(w, mus):
    w = np.asarray(w, complex); lg = np.zeros_like(w)
    for m in mus: lg += -(w + m) / 2 * LNPI + loggamma((w + m) / 2) + (1 - w + m) / 2 * LNPI - loggamma((1 - w + m) / 2)
    return np.exp(lg)

# local types: 1/L_l(X) coefficients and Bhargava masses (derivation 2, 3; re-derived: char poly of Frobenius on the permutation
# representation divided by (1 - X), ramified ones from the residue degrees)
TY = {'1111': ([1, -3, 3, -1], lambda l: 1 / 24), '211': ([1, -1, -1, 1], lambda l: 1 / 4), '22': ([1, 1, -1, -1], lambda l: 1 / 8),
      '31': ([1, 0, 0, -1], lambda l: 1 / 3), '4': ([1, 1, 1, 1], lambda l: 1 / 4), '1^2 11': ([1, -2, 1], lambda l: 1 / (2 * l)),
      '1^2 2': ([1, 0, -1], lambda l: 1 / (2 * l)), '1^2 1^2': ([1, -1], lambda l: 1 / (2 * l * l)), '2^2': ([1, 1], lambda l: 1 / (2 * l * l)),
      '1^3 1': ([1, -1], lambda l: 1 / l**2), '1^4': ([1], lambda l: 1 / l**3)}
def pol(c, X): return sum(ck * X**k for k, ck in enumerate(c))
def Ftype(nm, l, X): c = TY[nm][0]; return l / (l - 1) * pol(c, 1 / (l * X)) / pol(c, X)
def E_lim(l, X):
    Z = 1 + 1 / l + 2 / l**2 + 1 / l**3; return sum(TY[n][1](l) * Ftype(n, l, X) for n in TY) / Z
def E_obs(l, X, fr): return sum(f * Ftype(n, l, X) for n, f in fr.items())
def P_l(l, X):   # closed local factor / b_l(0): (1 - X^2/l) h_l(X)
    return (1 - X**2 / l) * (1 - 1 / (X * l**2) + X**2 / l - X**4 / l - X**5 / l + X**5 / l**2 - X**6 / l + X**7 / l**2)
# check the closed form against the type average at a few points: E_l = b_l(0) P_l / ((1-X^2)(1-X^3)(1-X^4))
for l in (2, 3, 11):
    for X in (0.2 + 0.1j, 0.37):
        b0 = 1 / (1 + l**-2 - l**-3 - l**-4); assert abs(E_lim(l, X) - b0 * P_l(l, X) / ((1 - X**2) * (1 - X**3) * (1 - X**4))) < 1e-12

PB = sieve(10**7).astype(float)
A1 = float(np.exp(-np.sum(np.log1p(PB**-2 - PB**-3 - PB**-4))))
def Qreal(x):   # prod_l P_l(l^-x), accelerated by the leading term -l^(x-2)
    X = PB**-x; return float(np.exp(np.sum(np.log(P_l(PB, X) / (1 - PB**(x - 2)))))) / float(szeta(2 - x))
res2 = A1 * szeta(1.5) * szeta(2.0) * Qreal(0.5) / 2; res3 = A1 * szeta(2 / 3) * Qreal(1 / 3) * szeta(4 / 3) / 3
print('A1 = %.6f; residue of D at 1/2 = %.6f (claimed 1/2); at 1/3 = %.5f' % (A1, res2, res3), flush=True)


t = np.arange(0, T + DT / 2, DT); w = C0 + 1j * t
PS = sieve(3000).astype(float); lr = np.zeros(len(w), complex)
for l in PS:
    X = np.exp(-w * math.log(l)); lr += np.log(P_l(l, X) / (1 - 1 / (X * l * l)))
Dw = A1 * zeta_c(2 * w) * zeta_c(3 * w) * zeta_c(4 * w) * np.exp(lr) / zeta_c(2 - w)
SIG = {0: (0, 0, 0), 1: (0, 0, 1), 2: (0, 1, 1)}
def M1T(mus, y, Tc):
    keep = t <= Tc + 1e-9; tt = t[keep]; ww = w[keep]; sw = np.ones(len(tt)); sw[1:-1:2] = 4; sw[2:-1:2] = 2; sw *= DT / 3
    g = Dw[keep] * Kt(ww, mus) / (1 - ww) * sw; out = np.zeros(len(y))
    for a in range(0, len(y), 40): out[a:a+40] = (np.exp(np.outer(np.log(y[a:a+40]), 1 - ww)) @ g).real / math.pi
    for w0, r in [(0.5, res2), (1/3, res3)]: out += r * float(Kt(np.array([w0 + 0j]), mus)[0].real) * y**(1 - w0) / (1 - w0)
    return out - np.sqrt(y)
for y0 in (1/27, 0.25, 1.0, 2.25):
    e = y0 + np.arange(-6, 7) * 0.004 * y0
    for r2 in (1, 2):
        rows = []
        for Tc in (500.0, 1000.0, 2000.0):
            F = M1T(SIG[r2], e, Tc); rows.append(np.diff(F) / np.diff(e))
        print('y0=%.4f r2=%d bins of width %.4f:' % (y0, r2, 0.004 * y0))
        for Tc, rr in zip((500, 1000, 2000), rows): print('   T=%4d %s' % (Tc, ' '.join('%6.2f' % v for v in rr)))
