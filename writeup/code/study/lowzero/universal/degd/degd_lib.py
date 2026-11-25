"""Library for the general-degree one-point identities (G1, G2, G3).
Analytic normalisation: L(s) = sum lam(n) n^-s, Lambda(s) = Q^{s/2} prod_j Gamma_R(s+mu_j) L(s) = eps * conj(Lambda(1 - conj s)).
X(s) = Q^{1/2-s} prod_j Gamma_R(1-s+conj mu_j)/Gamma_R(s+mu_j);  X(1/2+it)^{-1} = e^{2 i theta(t)} with
theta(t) = (t/2) ln Q + sum_j [Im lnGamma((1/2+mu_j+it)/2) - ((t + Im mu_j)/2) ln pi].
Window: H = 1_[T1,T2] * g_Delta, g_Delta(z) = exp(-z^2/Delta^2)/(sqrt(pi) Delta); Hhat(xi) = int H e^{-i xi t} dt.
"""
import math, cmath
import numpy as np
from scipy.special import loggamma, digamma, erf

PI = math.pi

def theta(t, Q, mu):
    t = np.asarray(t, float)
    out = 0.5 * t * math.log(Q)
    for m in mu:
        m = complex(m)
        out = out + np.imag(loggamma((0.5 + m + 1j * t) / 2)) - 0.5 * (t + m.imag) * math.log(PI)
    return out

def theta1(t, Q, mu):
    """theta'(t)."""
    t = np.asarray(t, float)
    out = 0.5 * math.log(Q) + 0 * t
    for m in mu:
        m = complex(m)
        out = out + 0.5 * np.real(digamma((0.5 + m + 1j * t) / 2)) - 0.5 * math.log(PI)
    return out

def H_box(t, T1, T2, D):
    t = np.asarray(t, float)
    return 0.5 * (erf((t - T1) / D) - erf((t - T2) / D))

def Hhat_box(xi, T1, T2, D):
    xi = np.asarray(xi, float)
    out = np.empty(xi.shape, complex)
    small = np.abs(xi) < 1e-12
    xs = np.where(small, 1.0, xi)
    W = (np.exp(-1j * xs * T1) - np.exp(-1j * xs * T2)) / (1j * xs)
    W = np.where(small, T2 - T1, W)
    return W * np.exp(-D * D * xi * xi / 4)

def grid(T1, T2, D, h=0.05, cut=6.2):
    a, b = T1 - cut * D, T2 + cut * D
    n = int(math.ceil((b - a) / h))
    return a + h * np.arange(n + 1), h

def osc_sum(lnu, tg, f, h, sign=-1, block=200000):
    """tg must be tg[0] + h*arange(K)."""
    """I(u) = h * sum_k f_k exp(sign*i*t_k*ln u) for every entry of lnu (trapezoid; f vanishes at both ends).
    Recurrence in k with an exact restart every 128 steps."""
    lnu = np.asarray(lnu, float)
    out = np.zeros(lnu.shape, complex)
    K = len(tg)
    for s0 in range(0, len(lnu), block):
        L = lnu[s0:s0 + block]
        z = np.exp(sign * 1j * h * L)
        acc = np.zeros(L.shape, complex)
        k = 0
        while k < K:
            v = np.exp(sign * 1j * tg[0] * L) * np.exp(sign * 1j * (h * k) * L)
            kend = min(K, k + 128)
            for j in range(k, kend):
                acc += f[j] * v
                v *= z
            k = kend
        out[s0:s0 + block] = h * acc
    return out

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def prime_powers(N):
    """arrays (n, p, k) of all prime powers n = p^k <= N, sorted by n."""
    P = primes_upto(N)
    ns, ps, ks = [], [], []
    for p in P:
        p = int(p); q = p; k = 1
        while q <= N:
            ns.append(q); ps.append(p); ks.append(k); q *= p; k += 1
    o = np.argsort(ns)
    return np.array(ns)[o], np.array(ps)[o], np.array(ks)[o]

def power_sums_from_poly(P, K):
    """P = [1, c1, ..., cd] (L_p = 1/P(X), P(X) = prod (1 - a_j X)); return s_k = sum a_j^k, k = 1..K (Newton)."""
    e = [1.0 + 0j]
    d = len(P) - 1
    # P(X) = sum c_j X^j with c_j = (-1)^j e_j
    ee = [complex(P[j]) * (-1) ** j for j in range(d + 1)]
    s = [0j] * (K + 1)
    for k in range(1, K + 1):
        acc = 0j
        for j in range(1, min(k - 1, d) + 1):
            acc += (-1) ** (j - 1) * ee[j] * s[k - j]
        if k <= d:
            acc += (-1) ** (k - 1) * k * ee[k]
        s[k] = acc
    return s[1:]

def series_inv(P, J):
    """coefficients of 1/P(X) up to X^J."""
    lam = [0j] * (J + 1); lam[0] = 1 / complex(P[0])
    for j in range(1, J + 1):
        s = 0j
        for i in range(1, min(j, len(P) - 1) + 1):
            s += complex(P[i]) * lam[j - i]
        lam[j] = -s / complex(P[0])
    return lam

def multiplicative_from_local(N, local):
    """lam(n), n <= N (index n), from local(p) -> list lam(p^k), k = 0..K (must cover p^k <= N)."""
    lam = np.ones(N + 1, complex); lam[0] = 0
    P = primes_upto(N)
    for p in P:
        p = int(p)
        loc = local(p)
        q = p; k = 1
        while q <= N:
            # multiply lam at multiples of q not of q*p by loc[k]/loc[k-1]... do it directly: indices with v_p = k
            idx = np.arange(q, N + 1, q)
            if q * p <= N:
                idx = idx[(idx // q) % p != 0]
            lam[idx] *= loc[k]
            q *= p; k += 1
    return lam
