"""Independent check of the family criterion (mpmath and sympy only; shares the zero list and the budget delta).

1. A few zeros of the zero list, refined again with mpmath's own Dirichlet L-function (Hurwitz sums at 20 digits).
2. The total number of zeros against the main term of the counting formula (sanity, not proof).
3. Lambda_hat(n) recomputed in mpmath for a few n: own primitive root and discrete logarithms (pow), own
   Kaiser-Bessel transform (mpmath Bessel functions), own P(x) and trivial-zero integrals (mp.quad).
Usage: fam_check.py zeros.npz famcrit.npy c beta m
"""
import sys, math
import numpy as np
import mpmath as mp
from sympy import primitive_root, isprime, factorint

mp.mp.dps = 20
zf, rf = sys.argv[1], sys.argv[2]; c, beta, m = mp.mpf(sys.argv[3]), mp.mpf(sys.argv[4]), int(sys.argv[5])
D = np.load(zf); q = int(D['q']); T = float(D['T']); kz = D['k']; gz = D['gamma']
R = np.load(rf)
g = int(primitive_root(q))
logtab = {}
v = 1
for j in range(q - 1):
    logtab[v] = j; v = v * g % q
out = []
def say(s):
    print(s, flush=True); out.append(s)

# 1
rng = np.random.default_rng(5)
for i in rng.choice(len(gz), 3, replace=False):
    k = int(kz[i]); gam = float(gz[i])
    w = mp.exp(2j * mp.pi / (q - 1))
    chi = [0] + [w ** ((logtab[a] * k) % (q - 1)) for a in range(1, q)]
    z = mp.findroot(lambda s: mp.dirichlet(s, chi), mp.mpc(0.5, gam))
    say(f'1. chi_{k}: zero {gam:.15f} in the list, mpmath {mp.nstr(z.imag, 18)} (real part {mp.nstr(z.real, 6)}); difference {float(abs(z.imag - gam)):.1e}')

# 2
cnt = int(D['count'][1:].sum())
main = (q - 2) * T / math.pi * math.log(q * T / (2 * math.pi * math.e)) + 0.25
say(f'2. zeros with |gamma| <= T over the {q - 2} L-functions: {cnt}; main term of the counting formula {main:.0f}; '
    f'difference per L-function {(cnt - main) / (q - 2):+.4f} (the bound allows {0.22737 * math.log(q * (T + 2) / (2 * math.pi)) + 2 * math.log(1 + math.log(q * (T + 2) / (2 * math.pi))) - 0.25:.2f})')

# 3
C0 = mp.sqrt(2 * mp.pi) * beta ** m
def Zmp(s):
    nu = m + mp.mpf(1) / 2
    if s < beta:
        qq = mp.sqrt(beta ** 2 - s ** 2); return mp.besseli(nu, qq) / qq ** nu
    y = mp.sqrt(s ** 2 - beta ** 2); return mp.besselj(nu, y) / y ** nu
W = [2 * Zmp(c * mp.mpf(float(t))) for t in gz]
k0 = mp.besseli(m, beta)
def Kmp(v, a):
    t = 1 - (v / a) ** 2
    return t ** (mp.mpf(m) / 2) * mp.besseli(m, beta * mp.sqrt(t)) if t > 0 else mp.mpf(0)
for n in (3, 97, 1009, 4096, 9973, 9999):
    x = mp.mpf(n); a = c * x; ln = mp.log(x); idn = logtab[n % q]
    S = mp.fsum(Wz * mp.cos(mp.mpf(float(t)) * ln - 2 * mp.pi * ((idn * int(k)) % (q - 1)) / (q - 1)) for Wz, t, k in zip(W, gz, kz))
    nu = m + mp.mpf(1) / 2; r = a / (2 * x); qq = mp.sqrt(beta ** 2 + r ** 2)
    Px = a * C0 * mp.besseli(nu, qq) / qq ** nu
    tt = lambda v: x * mp.e ** (v / x)
    Ez = mp.sqrt(x) / x * mp.quad(lambda v: Kmp(v, a) * mp.sqrt(tt(v)) / (tt(v) ** 3 - tt(v)), [-a, 0, a])
    Ee = mp.sqrt(x) / x * mp.quad(lambda v: Kmp(v, a) * tt(v) ** 1.5 / (tt(v) ** 2 - 1), [-a, 0, a])
    lam_hat = (Px - Ez + Ee - c * mp.sqrt(x) * C0 * S) / ((q - 1) * k0)
    f = factorint(n); lam = math.log(next(iter(f))) if len(f) == 1 else 0.0
    row = R[R[:, 0] == n][0]
    say(f'3. n = {n}: Lambda_hat in mpmath {mp.nstr(lam_hat, 12)}, in the run {row[1]:.12f}, Lambda(n) = {lam:.6f}, delta = {row[3]:.4f}; '
        f'inside: {abs(lam_hat - lam) <= row[3]}; prime: {isprime(n)}, decided prime: {bool(lam_hat > 0.75 * math.log(n))}')
open(rf.replace('.npy', '_check.txt'), 'w').write('\n'.join(out) + '\n')
