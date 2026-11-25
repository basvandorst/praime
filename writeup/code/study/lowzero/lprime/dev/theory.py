"""Theory numbers for log|L'(rho, chi_D)| (derivation.md). Prints:
 1. CUE(N): exact cumulants and moments of log|Lambda'(theta_n)| (Palm measure = Keating-Snaith moment with s = 1 + k);
 2. one Euler factor p (uniform phase): cumulants kappa_1..4 of A_p = -ln|1 - p^(-1/2) e^(i phi)|, and the first-order
    (1/L) correction of the local-density-modulation (LDM) model, per unit ln p / L;
 3. bad-prime moment factors r_p(k) = 1/2F1(k,k;1;1/p) and their 1/L corrections (2 ln p/L)(1+k)^2 G_p(2k);
 4. the left-tail constant of CUE and its bad-prime factor.
"""
import numpy as np, mpmath as mp
from math import log, pi
mp.mp.dps = 30
EG = float(mp.euler)

# ---------- 1. CUE ----------
def lnM(n, s):  # ln E_CUE(n) |Lambda(1)|^(2s), Keating-Snaith, analytic in n and s
    return (mp.log(mp.barnesg(n + 1)) + mp.log(mp.barnesg(n + 2 * s + 1)) + 2 * mp.log(mp.barnesg(1 + s))
            - 2 * mp.log(mp.barnesg(n + s + 1)) - mp.log(mp.barnesg(1 + 2 * s)))

def cue_K(N, t):  # cumulant generating function of Y = log|Lambda'(theta_n)| in CUE(N): E e^{tY} = M_{N-1}(1+t/2)/M_{N-1}(1)
    return lnM(N - 1, 1 + t / 2) - lnM(N - 1, 1)

def cue_cumulants(N, nmax=4):
    return [float(mp.diff(lambda t: cue_K(N, t), 0, j)) for j in range(1, nmax + 1)]

def cue_moment(N, k):  # E |Lambda'|^{2k}
    return float(mp.exp(cue_K(N, 2 * k)))

# ---------- 2. Euler factor ----------
M = 1 << 14; PHI = 2 * pi * (np.arange(M) + 0.5) / M
def series_exp_moments(f, w, order):  # coefficients c_j of E[w e^{t f}] = sum c_j t^j
    out = []; fj = np.ones_like(f); fact = 1.0
    for j in range(order + 1):
        out.append(float(np.mean(w * fj)) / fact); fj = fj * f; fact *= (j + 1)
    return np.array(out)

def ser_mul(a, b, n): return np.array([sum(a[i] * b[j - i] for i in range(j + 1)) for j in range(n + 1)])
def ser_inv(a, n):
    b = np.zeros(n + 1); b[0] = 1 / a[0]
    for j in range(1, n + 1): b[j] = -sum(a[i] * b[j - i] for i in range(1, j + 1)) / a[0]
    return b
def ser_log(a, n):  # log of series with a[0] = 1
    d = np.array([(j + 1) * a[j + 1] for j in range(n)]); q = ser_mul(d, ser_inv(a, n), n - 1)
    return np.concatenate([[0.0], [q[j - 1] / j for j in range(1, n + 1)]])

def euler_terms(p, order=6):
    """Return (lead, corr): series in t of K0(t) = ln E e^{tA} and of K1(t)/(ln p/L) =
    -2 (1+t/2)^2 E[e^{tA} R]/E[e^{tA}], R = Re z/(1-z), z = p^{-1/2} e^{i phi}."""
    z = p ** -0.5 * np.exp(1j * PHI); A = -np.log(np.abs(1 - z)); R = np.real(z / (1 - z))
    one = np.ones(M); E0 = series_exp_moments(A, one, order); ER = series_exp_moments(A, R, order)
    K0 = ser_log(E0, order); G = ser_mul(ER, ser_inv(E0, order), order)
    pre = np.zeros(order + 1); pre[0] = 1; pre[1] = 1; pre[2] = 0.25
    K1 = -2 * ser_mul(pre, G, order)
    fact = np.array([float(mp.factorial(j)) for j in range(order + 1)])
    return K0 * fact, K1 * fact  # cumulants kappa_j = j! * coefficient

def hyp(k, p): return float(mp.hyp2f1(k, k, 1, mp.mpf(1) / p))

def G_of(p, t):  # G(t) = E[|1-z|^{-t} R] / E[|1-z|^{-t}]
    z = p ** -0.5 * np.exp(1j * PHI); w = np.abs(1 - z) ** (-t); R = np.real(z / (1 - z))
    return float(np.mean(w * R) / np.mean(w))

if __name__ == '__main__':
    print('== 1. CUE(N): cumulants of log|Lambda\'(theta_n)| (exact) ==')
    print('  N      mean   (H_N-1)   var      k3        k4      E|.|^-2  E|.|^-1   E|.|^1    E|.|^2')
    for N in (5, 8, 10, 12, 16, 20, 40, 100):
        c = cue_cumulants(N); HN = float(mp.harmonic(N)) - 1
        print('%4d  %8.5f %8.5f %8.5f %9.5f %9.5f  %8.5f %8.5f %8.4f %9.3f' % (N, c[0], HN, c[1], c[2], c[3], cue_moment(N, -1), cue_moment(N, -0.5), cue_moment(N, 0.5), cue_moment(N, 1)))
    print('  large-N: mean = ln N + gamma - 1 + 1/(2N) + ...; var = 1/2 ln N + (gamma+1)/2 + 1 - pi^2/4 = 1/2 ln N %+.5f' % ((EG + 1) / 2 + 1 - pi ** 2 / 4))
    print('  k3(inf), k4(inf) = ', [round(float(mp.diff(lambda t: cue_K(400, t), 0, j)), 5) for j in (3, 4)])
    print()
    print('== 2. one Euler factor (good prime): cumulants of A_p, and 1/L corrections (per unit ln p/L) ==')
    print('  p   k1   var=Li2/2   k3        k4       |  dk1       dvar       dk3       dk4    (multiply by ln p/L)')
    for p in (2, 3, 5, 7, 11, 13):
        K0, K1 = euler_terms(p)
        print('%3d  %.1e %9.5f %9.5f %9.5f | %9.5f %9.5f %9.5f %9.5f   Li2/2 check %.5f, ln(p/(p-1)) %.5f' % (p, K0[1], K0[2], K0[3], K0[4], K1[1], K1[2], K1[3], K1[4], float(mp.polylog(2, mp.mpf(1) / p)) / 2, log(p / (p - 1))))
    print()
    print('== 3. bad-prime factor on E|L\'|^{2k}: r_p(k) [1 + (2 ln p/L)(1+k)^2 G_p(2k)] ==')
    for p in (2, 3, 5, 7, 11, 13):
        row = []
        for k in (-1.5, -1, -0.5, 0.5, 1):
            row.append('k=%+.1f r=%.5f c=%+.4f' % (k, 1 / hyp(k, p), 2 * log(p) * (1 + k) ** 2 * G_of(p, 2 * k)))
        print('%3d ' % p + ' | '.join(row))
    print('  checks: r_p(-1) = p/(p+1): %s; r_p(1) = 1-1/p: %s; G_p(2) = 1/(p-1): %s' % (
        [round(1 / hyp(-1, p) - p / (p + 1), 12) for p in (2, 3, 5)], [round(1 / hyp(1, p) - (1 - 1 / p), 12) for p in (2, 3, 5)],
        [round(G_of(p, 2) - 1 / (p - 1), 10) for p in (2, 3, 5)]))
