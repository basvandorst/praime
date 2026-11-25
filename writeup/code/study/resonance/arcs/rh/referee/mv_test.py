"""Referee test of the hybrid mean value (Lemma 5 of rh_proof.md) and of the published forms.

LHS(c) = sum_{|m|<=M} int_{-T}^{T} |sum_a c(a) lambda^m(a) Na^{-it}|^2 dt = c^* K c,
K_ab = D_M(4(th_a - th_b)) * K_T(log(Na/Nb)),
D_M(x) = sum_{|m|<=M} e^{imx},   K_T(u) = int_{-T}^{T} e^{iut} dt = 2 sin(Tu)/u.
lambda(a) = (alpha/|alpha|)^4, so lambda^m(a) = exp(4 i m th_a).

Bounds compared:
  L5   : e^2 sqrt(pi) (sqrt(pi) M + pi N/2) sum |c|^2 (T + sqrt(pi) D(a)),  D(a) = cont(a) sqrt(N/Na)
  ST   : (M T + N min(M, T)) sum |c|^2        (Stucky Lemma 4.2, constant 1 shown; the claim is <<)
  DUKE : sum (Na + M T) |c|^2                 (Duke 1989 p. 205 conjecture, q = 1, n = 2)
Also the exact worst case over ALL coefficient vectors: lambda_max of W^{-1/2} K W^{-1/2}.
"""
import sys
import math
import time
import numpy as np
import scipy.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ideals import enumerate_ideals

E2SP = math.e ** 2 * math.sqrt(math.pi)


def kernel(I, idx, M, T):
    th = I['theta'][idx]
    ln = np.log(I['norm'][idx].astype(float))
    pr = I['prim'][idx]
    dth = th[:, None] - th[None, :]
    y = np.mod(dth + math.pi / 4, math.pi / 2) - math.pi / 4  # in [-pi/4, pi/4)
    same = pr[:, None] == pr[None, :]
    s = np.sin(2 * y)
    with np.errstate(divide='ignore', invalid='ignore'):
        D = np.sin((2 * M + 1) * 2 * y) / s
    D[same] = 2 * M + 1
    # non-collinear pairs must have |sin| bounded below (spacing check)
    if (~same).any():
        assert np.abs(s[~same]).min() > 0.5 / (I['norm'][idx].max() ** 1.0) * 0.9
    u = ln[:, None] - ln[None, :]
    KT = 2 * T * np.sinc(T * u / math.pi)
    return D * KT


def lhs(K, c):
    return float(np.real(np.conj(c) @ (K @ c)))


def bounds(I, idx, c, M, T, N):
    n = I['norm'][idx].astype(float)
    Dw = I['cont'][idx] * np.sqrt(N / n)
    c2 = np.abs(c) ** 2
    L5 = E2SP * (math.sqrt(math.pi) * M + math.pi * N / 2) * float(np.sum(c2 * (T + math.sqrt(math.pi) * Dw)))
    ST = (M * T + N * min(M, T)) * float(np.sum(c2))
    DU = float(np.sum((n + M * T) * c2))
    return L5, ST, DU


def brute_check():
    """Direct evaluation (sum over m, Gauss-Legendre in t) against the kernel form."""
    I = enumerate_ideals(60)
    rng = np.random.default_rng(1)
    idx = np.arange(len(I['norm']))
    c = rng.normal(size=len(idx)) + 1j * rng.normal(size=len(idx))
    M, T = 7, 3.0
    K = kernel(I, idx, M, T)
    v1 = lhs(K, c)
    x, w = np.polynomial.legendre.leggauss(400)
    t = T * x
    w = T * w
    tot = 0.0
    lam = np.exp(4j * I['theta'])
    for m in range(-M, M + 1):
        coef = c * lam ** m
        P = (coef[None, :] * np.exp(-1j * np.outer(t, np.log(I['norm'].astype(float))))).sum(1)
        tot += float(np.sum(w * np.abs(P) ** 2))
    return v1, tot


def coeff_families(I, idx, rng):
    n = I['norm'][idx]
    cont = I['cont'][idx]
    a = I['a'][idx]
    b = I['b'][idx]
    fam = {}
    fam['rational (k)'] = (b == 0).astype(complex)
    big = cont >= max(2, int(math.isqrt(int(n.max())) // 6))
    fam['content >= sqrt(N)/6'] = big.astype(complex)
    # one ray: multiples d(2+i)
    ray = (I['prim'][idx] == 2 * 100003 + 1)
    fam['ray of 2+i'] = ray.astype(complex)
    # thin angular window near a generic angle, all contents
    th = I['theta'][idx]
    th0 = 0.61803
    win = np.abs(th - th0) < 2.0 / n.max() ** 0.5
    fam['angle window 2/sqrt(N)'] = win.astype(complex)
    fam['primes'] = I['isprime'][idx].astype(complex)
    fam['random signs'] = rng.choice([-1.0, 1.0], size=len(idx)).astype(complex)
    fam['random phases'] = np.exp(2j * math.pi * rng.random(len(idx)))
    fam['divisor d(a)'] = I['d'][idx].astype(complex)
    fam['Moebius'] = I['mu'][idx].astype(complex)
    fam['ones'] = np.ones(len(idx), dtype=complex)
    fam['log N'] = np.log(n.astype(float)).astype(complex)
    return fam


def main():
    t0 = time.time()
    v1, v2 = brute_check()
    print(f"brute-force check of the kernel form: kernel {v1:.10e}  direct {v2:.10e}  rel {abs(v1 - v2) / abs(v2):.2e}")
    assert abs(v1 - v2) / abs(v2) < 1e-8
    rng = np.random.default_rng(7)
    regimes = [(1, 1.0), (10, 1.0), (300, 1.0), (3000, 1.0), (100000, 1.0), (100000, 0.1),
               (1, 1000.0), (30, 30.0), (1000, 1000.0), (10, 10000.0), (1000, 10.0), (3, 0.3)]
    worst_family = 0.0
    worst_eig = 0.0
    for N in (500, 2500):
        I = enumerate_ideals(N)
        idx = np.arange(len(I['norm']))
        fams = coeff_families(I, idx, rng)
        print(f"\nN = {N}: {len(idx)} ideals with norm <= N")
        for (M, T) in regimes:
            K = kernel(I, idx, M, T)
            # worst case over all coefficients for Lemma 5
            n = I['norm'][idx].astype(float)
            w = E2SP * (math.sqrt(math.pi) * M + math.pi * N / 2) * (T + math.sqrt(math.pi) * I['cont'][idx] * np.sqrt(N / n))
            s = 1 / np.sqrt(w)
            Kw = K * s[:, None] * s[None, :]
            lam = sla.eigh(Kw, eigvals_only=True, subset_by_index=[len(idx) - 1, len(idx) - 1])[0]
            lamK = sla.eigh(K, eigvals_only=True, subset_by_index=[len(idx) - 1, len(idx) - 1])[0]
            worst_eig = max(worst_eig, lam)
            line = []
            for name, c in fams.items():
                if not np.any(c):
                    continue
                L = lhs(K, c)
                L5, ST, DU = bounds(I, idx, c, M, T, N)
                worst_family = max(worst_family, L / L5)
                line.append((name, L / L5, L / ST, L / DU))
            print(f"  M={M:>6} T={T:>7}: max over ALL c of LHS/L5 = {lam:.4f};  max LHS/(sum|c|^2) = {lamK:.4e}"
                  f" vs MT+N min(M,T) = {M*T + N*min(M, T):.4e} (ratio {lamK/(M*T + N*min(M, T)):.2f})")
            for name, r5, rs, rd in line:
                print(f"      {name:<24} LHS/L5 {r5:8.4f}   LHS/ST {rs:9.3f}   LHS/DUKE {rd:9.3f}")
    print(f"\nlargest LHS/L5 over the families: {worst_family:.4f}; over all c (eigen): {worst_eig:.4f}")
    assert worst_eig <= 1.0
    print(f"time {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()
