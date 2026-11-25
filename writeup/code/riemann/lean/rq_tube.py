# The prime sum of the one-sided criterion (grh-lean, New/Reduction/Zeros/Criterion.lean), on primes, for q = 8, d = 7:
#   R_q(x) = sum_{n <= x, n = 1 mod q} Lambda(n) (1 - (n/d)) (1/n - 1/x) - log x / phi(q) - c_q,
#   c_q = Re G_q(1), G_q = (sum_psi (L'/L(psi chi_d) - L'/L(psi)) - 1)/phi(q) with the L-functions made entire.
# Explicit formula (residues of the transform: (n(s) - n'(s))/(phi s (1 - s))):
#   R_q(x) = (1/phi) [ sum over zeros of L(psi) - sum over zeros of L(psi chi_d) ] x^(rho-1)/(rho (1 - rho)) + O(1/x).
# With GRH for Dirichlet L-functions: |R_q(x)| sqrt x <= C = (1/phi) sum over all these zeros of 1/(rho (1 - rho)).
import math, sys, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(400_000_000)
q, d, phi = 8, 7, 4
PLAIN = [1, -4, 8, -8]          # primitive characters inducing the characters mod 8 (Kronecker symbols)
TWIST = [-7, 28, -56, 56]       # the same times (./7) = Kronecker symbol of -7
g = 0.5772156649015329
def ldl(D):                      # L'/L(1, chi_D), D != 1
    return float(pari('lfun(%d, 1, 1)/lfun(%d, 1)' % (D, D)))
def zsum(D):                     # sum over zeros of 1/(rho (1 - rho)) = log(|D|/pi) + psi((1 + a)/2) + 2 L'/L(1, chi_D)
    if D == 1: return 2 + g - math.log(4 * math.pi)
    a = 1 if D < 0 else 0
    return math.log(abs(D) / math.pi) + float(pari('psi(%s)' % ('1' if a else '1/2'))) + 2 * ldl(D)
# entire L of the characters mod 8 and mod 56: Euler factors lost at 2 (for D = 1 and D = -7) add log 2 at s = 1
G1 = ((ldl(-7) + math.log(2)) + ldl(28) + ldl(-56) + ldl(56) - (g + math.log(2)) - ldl(-4) - ldl(8) - ldl(-8) - 1) / phi
C = sum(zsum(D) for D in PLAIN + TWIST) / phi
print('c_q = Re G_q(1) = %.10f   tube constant C = %.10f' % (G1, C))
# check zsum against zeros for one character
z = [float(t) for t in pari('lfunzeros(lfuncreate(-7), 400)')]
print('check: chi_-7, zeros to 400: 2 sum 1/(1/4 + g^2) = %.6f, closed form %.6f (tail about 2 int_400^inf log(7 t/2pi)/(2 pi t^2) dt = %.6f)' % (
    2 * sum(1 / (0.25 + t * t) for t in z), zsum(-7), (math.log(7 * 400 / (2 * math.pi)) + 1) / (math.pi * 400)))
# the prime sum up to X
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**8
SEG = 10**7
base = np.ones(int(math.isqrt(X)) + 2, bool); base[:2] = False
for i in range(2, int(math.isqrt(len(base))) + 1):
    if base[i]: base[i*i::i] = False
small = np.nonzero(base)[0]
qr = np.zeros(7, bool); qr[[1, 2, 4]] = True
ns, ws = [], []
for lo in range(0, X + 1, SEG):
    hi = min(lo + SEG, X + 1); s = np.ones(hi - lo, bool)
    if lo == 0: s[:2] = False
    for p in small:
        if p * p >= hi: break
        st = max(p * p, (lo + p - 1) // p * p); s[st - lo::p] = False
    pr = np.nonzero(s)[0] + lo
    m = (pr % 8 == 1) & (pr != 7)
    pp = pr[m]; lam = np.log(pp) * (1.0 - np.where(qr[pp % 7], 1.0, -1.0))
    ns.append(pp[lam > 0]); ws.append(lam[lam > 0])
for p in small:                  # prime powers p^k, k >= 2, = 1 mod 8
    p = int(p); pk = p * p; k = 2
    while pk <= X:
        if pk % 8 == 1:
            chi = 0.0 if p == 7 else (1.0 if qr[pk % 7] else -1.0)         # (n/7) = 0 on powers of 7: weight log 7
            if chi < 1: ns.append(np.array([pk])); ws.append(np.array([(1 - chi) * math.log(p)]))
        pk *= p; k += 1
n = np.concatenate(ns); w = np.concatenate(ws); o = np.argsort(n); n = n[o].astype(float); w = w[o]
A = np.cumsum(w / n); B = np.cumsum(w)
xs = np.exp(np.linspace(math.log(100), math.log(X), 4001))
idx = np.searchsorted(n, xs, side='right') - 1
R = A[idx] - B[idx] / xs - np.log(xs) / phi - G1
T = R * np.sqrt(xs)
print('terms: %d prime powers = 1 mod 8 with (n/7) != 1, up to %d' % (len(n), X))
for a_, b_ in ((1e2, 1e4), (1e4, 1e6), (1e6, X)):
    k = (xs >= a_) & (xs <= b_)
    print('  x in [%.0e, %.0e]: R sqrt(x) between %+.4f and %+.4f' % (a_, b_, T[k].min(), T[k].max()))
print('largest |R| sqrt(x) on the grid: %.4f = %.3f of C ; R(X) sqrt(X) = %+.4f' % (np.abs(T).max(), np.abs(T).max() / C, T[-1]))
# explicit formula from zeros at a few x
zs = {D: np.array([float(t) for t in pari('lfunzeros(lfuncreate(%d), 600)' % D)]) for D in PLAIN + TWIST}
for x in (1e4, 1e6, 1e8):
    if x > X: continue
    tot = 0.0
    for D in PLAIN + TWIST:
        sg = 1 if D in PLAIN else -1
        tot += sg * float(np.sum(2 * np.cos(zs[D] * math.log(x)) / (0.25 + zs[D] ** 2)))      # pairs rho, conj rho
    i = np.searchsorted(n, x, side='right') - 1
    Rx = (A[i] - B[i] / x - math.log(x) / phi - G1) * math.sqrt(x)
    print('  x = %.0e: R sqrt x from primes %+.5f ; from the zeros below 600: %+.5f' % (x, Rx, tot / phi))
