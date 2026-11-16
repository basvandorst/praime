"""Theory 11 (no sieve in the theory): the race between non-residue and residue primes mod a prime P is decided by the zeros of
L(s, chi_P). Explicit formula (GRH): E(x) = (ln x / sqrt x) (pi_N(x) - pi_R(x)) = 1 + sum_rho x^(i gamma)/rho; the 1 is the prime
squares. Over all zeros the terms have total variance V = sum_rho 1/|rho|^2 = 2 [Re L'/L(1, chi) + 1/2 ln(P/pi) + 1/2 psi((1+a)/2)]
(Hadamard product; a = 0 for even, 1 for odd chi). One-zero law: the lowest zero gamma_1 gives a wave of amplitude A = 2/|rho_1|,
the rest is noise of variance V - A^2/2, so the share of (logarithmic) time with E < 0 is P(1 + A cos(gamma_1 ln x + phase) + G < 0).
Predictions from L-function values only; then measured on the primes to 10^8."""
import numpy as np, mpmath as mp, math
mp.mp.dps = 20
PS = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
def chi_list(P):
    sq = {(i * i) % P for i in range(1, P)}
    return [0] + [1 if r in sq else -1 for r in range(1, P)]
def Lfun(s, chi):
    P = len(chi)
    return mp.mpf(P) ** (-s) * mp.fsum(chi[r] * mp.zeta(s, mp.mpf(r) / P) for r in range(1, P) if chi[r])
def Lam(t, P, chi, a):
    s = mp.mpc(0.5, t)
    return (mp.mpf(P) / mp.pi) ** ((s + a) / 2) * mp.gamma((s + a) / 2) * Lfun(s, chi)
res = {}
for P in PS:
    chi = chi_list(P); a = 0 if chi[P - 1] == 1 else 1
    g1 = None; t = 0.02; v = float(mp.re(Lam(t, P, chi, a)))
    while g1 is None and t < 16:
        t2 = t + 0.04; v2 = float(mp.re(Lam(t2, P, chi, a)))
        if v * v2 <= 0:
            g1 = float(mp.findroot(lambda u: mp.re(Lam(u, P, chi, a)), (t, t2), solver='bisect'))
        t, v = t2, v2
    imag_ok = abs(float(mp.im(Lam(3.3, P, chi, a)))) < 1e-8 * max(1, abs(float(mp.re(Lam(3.3, P, chi, a)))))
    with mp.workdps(40):
        h = mp.mpf(10) ** -8
        Lp, Lm = Lfun(1 + h, chi), Lfun(1 - h, chi)
        L1 = (Lp + Lm) / 2; dL1 = (Lp - Lm) / (2 * h)
    V = 2 * float(mp.re(dL1 / L1) + 0.5 * mp.log(mp.mpf(P) / mp.pi) + 0.5 * mp.digamma((1 + a) / mp.mpf(2)))
    A = 2 / math.sqrt(0.25 + g1 ** 2)
    res[P] = (g1, A, V, a, imag_ok)
    print(f'P {P:>2} ({"odd " if a else "even"}): lowest zero {g1:.4f}, wave amplitude {A:.3f}, total variance {V:.3f}, rest {V - A * A / 2:.3f}, Lambda real on the line: {imag_ok}', flush=True)
X = 100_000_000
s = np.ones(X // 2, bool); s[0] = False
for i in range(3, int(X ** 0.5) + 1, 2):
    if s[i // 2]: s[i * i // 2::i] = False
pr = 2 * np.nonzero(s)[0] + 1
x0 = 1000
print(f'\nshare of log-time in [{x0:,}, {X:,}] with more residue primes than non-residue primes (E < 0)')
print('  P   measured   Gauss(all zeros)   one-zero law   one-zero with its phase   the same two with the finite-x square term b(x)')
lxg = np.linspace(math.log(x0), math.log(X), 4000)
bx = np.array([float(mp.li(mp.e ** (l / 2))) for l in lxg]) * 0.5 * lxg / np.exp(lxg / 2)
print(f'   b(x) = (ln x / sqrt x) * Li(sqrt x) / 2 runs from {bx[0]:.3f} at {x0:,} to {bx[-1]:.3f} at {X:,} (limit 1)')
from math import erf
Phi = lambda z: 0.5 * (1 + erf(z / math.sqrt(2)))
out = []
for P in PS:
    g1, A, V, a, _ = res[P]
    chi = np.array(chi_list(P))
    q = pr[pr != P]
    race = np.cumsum(-chi[q % P])
    w = np.diff(np.log(np.append(q, X)))
    m = q >= x0
    meas = np.sum(w[m] * (race[m] < 0)) / np.sum(w[m])
    sd = math.sqrt(max(V - A * A / 2, 1e-9))
    th = np.linspace(0, 2 * math.pi, 2001)[:-1]
    one = np.mean([Phi((-1 - A * math.cos(t)) / sd) for t in th])
    gauss = Phi(-1 / math.sqrt(V))
    lx = np.linspace(math.log(x0), math.log(X), 4000)
    rho = complex(0.5, g1)
    wave = np.array([2 * (np.exp(1j * g1 * l) / rho).real for l in lx])
    win = np.mean([Phi((-1 - v) / sd) for v in wave])
    gb = np.mean([Phi(-b / math.sqrt(V)) for b in bx])
    wb = np.mean([Phi((-b - v) / sd) for b, v in zip(bx, wave)])
    out.append((P, meas, gauss, one, win, gb, wb))
    print(f'{P:>3} {meas:>10.3f} {gauss:>18.3f} {one:>14.3f} {win:>20.3f} {gb:>14.3f} {wb:>10.3f}', flush=True)
o = np.array(out)
for j, nm in ((2, 'Gauss'), (3, 'one-zero'), (4, 'one-zero with phase'), (5, 'Gauss with b(x)'), (6, 'one-zero with phase and b(x)')):
    print(f'{nm}: correlation with measured {np.corrcoef(o[:, 1], o[:, j])[0, 1]:.3f}, mean absolute error {np.mean(np.abs(o[:, 1] - o[:, j])):.3f}')
