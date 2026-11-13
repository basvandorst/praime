"""The two extremes of the family mod 395,137: the zero closest to s = 1/2, and the character whose first zero is highest.
Z_k(t) evaluated directly (Euler-Maclaurin Hurwitz values, one character) at points near the zero; the zero refined by
bisection; |L(1/2, chi)| and the character's values at small primes reported."""
import math, numpy as np
from scipy.special import loggamma, bernoulli
D = np.load('lfam395137f.npz'); q = int(D['q']); pw = D['pw']; kz = D['k']; gz = D['gamma']
N = q - 1; ind = np.empty(q, dtype=np.int64); ind[pw] = np.arange(N)
NT_, M_ = 30, 9; alpha = pw / q
LOG = np.log(np.arange(NT_)[None, :] + alpha[:, None]); Na = NT_ + alpha; LNa = np.log(Na)
B = bernoulli(2 * M_); coef = [float(B[2 * j]) / math.factorial(2 * j) for j in range(1, M_ + 1)]
def hurwitz_all(s):
    main = np.exp(-s * LOG).sum(axis=1); Nas = np.exp(-s * LNa)
    tail = Nas * Na / (s - 1) + 0.5 * Nas; rising = s; powr = Nas / Na
    for j in range(1, M_ + 1):
        tail = tail + coef[j - 1] * rising * powr; rising = rising * (s + 2 * j - 1) * (s + 2 * j); powr = powr / (Na * Na)
    return main + tail
jj = np.arange(N)
def make(k):
    ch = np.exp(2j * math.pi * ((k * jj) % N) / N)
    tau = ch @ np.exp(2j * math.pi * pw / q); kap = k % 2
    eps = tau / ((1j) ** kap * math.sqrt(q)); rot = np.exp(-0.5j * np.angle(eps))
    def L(s): return np.exp(-s * math.log(q)) * (ch @ hurwitz_all(s))
    def Z(t):
        th = 0.5 * t * math.log(q / math.pi) + np.imag(loggamma((0.5 + kap + 1j * t) / 2))
        return float(np.real(rot * np.exp(1j * th) * L(0.5 + 1j * t)))
    return L, Z
primes = [p for p in range(2, 60) if all(p % r for r in range(2, int(p ** 0.5) + 1))]
for k in (590, 25256):
    kb = (-k) % N
    own = np.sort(gz[kz == k])[:3]; conj = np.sort(gz[kz == kb])[:3]
    L, Z = make(k)
    print(f'chi_{k} (order {N // math.gcd(k, N)}, {"odd" if k % 2 else "even"}): lowest zeros of chi on t > 0 {own.round(7).tolist()}, of its conjugate {conj.round(7).tolist()}')
    print(f'  L(1/2, chi) = {L(0.5 + 0j):.6e},  |L(1/2, chi)| = {abs(L(0.5 + 0j)):.3e}')
    lo_t = min(own[0], conj[0])
    if k == 590:
        # the zero nearest the centre: which side, and refine
        sgn = 1 if own[0] < conj[0] else -1
        Zs = make(k if sgn == 1 else kb)[1]
        a_, b_ = 0.0, 2 * lo_t
        za = Zs(a_)
        for _ in range(50):
            m_ = 0.5 * (a_ + b_); zm = Zs(m_)
            if zm * za > 0: a_, za = m_, zm
            else: b_ = m_
        print(f'  zero of {"chi" if sgn == 1 else "its conjugate"} at t = {0.5 * (a_ + b_):.12e}; Z(0) = {Zs(0.0):.3e}, Z(2 t0) = {Zs(2 * lo_t):.3e}')
    re = [math.cos(2 * math.pi * ((k * ind[p]) % N) / N) for p in primes]
    print('  Re chi(p) for p < 60: ' + ', '.join(f'{p}:{r:+.2f}' for p, r in zip(primes, re)) + f';  mean {np.mean(re):+.3f}')
