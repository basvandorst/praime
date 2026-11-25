#!/usr/bin/env python3
"""mock_sum.py -- the prime-sum part of the linear statistic of low zeros, for ALL characters mod a prime q at once.

For a prime q, every non-principal character chi_k (k = 1..q-2; chi_k(g) = e(k/(q-1)), g the least primitive root)
and every support s in the grid, with the Fejer test function fhat_s(u) = (1 - |u|/s)_+:

    P_s(chi) = -(2/log q) * sum_{n <= q^s, (n,q)=1} Lambda(n) n^(-1/2) fhat_s(log n/log q) * Re chi(n)

(= W_f(chi) - Wbar_{a(chi)}, the zero sum minus its gamma-factor part; derivation.md section 1).
Method: bin the weights by the discrete log of n mod q (array B_s of length q-1), one FFT of length q-1 per s:
    S_s(k) = sum_l B_s[l] e(k l/(q-1)),   P_s(chi_k) = -(2/log q) Re S_s(k).
Primes from a segmented odd-only sieve (no array above ~10 MB besides B); prime powers p^k (k >= 2) added exactly.
Also computed in the same pass: the independent-phase ("diagonal", Steinhaus) cumulants kappa_2..4 prime by prime.

Output (npz): power sums sum_chi P^j, j = 0..8, over the family (all / even / odd), the diagonal cumulants, and
(--save-p) the full array P_s(chi_k).  Resumable: a checkpoint is written every --ckpt-every segments.

usage: mock_sum.py --q 10007 --smax 2.0 --out q10007.npz [--save-p] [--seg 10000000]
Only numpy is needed.
"""
import argparse, math, os, sys, time
import numpy as np


# ---------- small number theory (self-contained) ----------
def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_factors(n):
    f, d = [], 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        f.append(n)
    return f


def primitive_root(q):
    fs = prime_factors(q - 1)
    g = 2
    while any(pow(g, (q - 1) // r, q) == 1 for r in fs):
        g += 1
    return g


def dlog_table(q, g):
    """L[r] = discrete log of r base g, r = 1..q-1 (L[0] unused = -1)."""
    L = np.full(q, -1, dtype=np.int64)
    x = 1
    for j in range(q - 1):
        L[x] = j
        x = x * g % q
    return L


def simple_sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    s[4::2] = False
    for i in range(3, int(math.isqrt(n)) + 1, 2):
        if s[i]:
            s[i * i::2 * i] = False
    return np.nonzero(s)[0].astype(np.int64)


def odd_primes_in(lo, hi, base):
    """odd primes in [lo, hi) using base primes (all primes <= sqrt(hi))."""
    if lo % 2 == 0:
        lo += 1
    if hi <= lo:
        return np.zeros(0, dtype=np.int64)
    n = (hi - lo + 1) // 2
    seg = np.ones(n, dtype=bool)
    for p in base:
        if p == 2:
            continue
        pp = int(p) * int(p)
        if pp >= hi:
            break
        start = max(pp, ((lo + p - 1) // p) * p)
        if start % 2 == 0:
            start += p
        if start >= hi:
            continue
        seg[(start - lo) // 2::int(p)] = False
    nums = lo + 2 * np.nonzero(seg)[0].astype(np.int64)
    return nums[nums > 2] if lo <= 1 else nums  # removes 1


def fhat_tri(u, s):
    return np.clip(1.0 - np.abs(u) / s, 0.0, None)


def diag_cumulants_one_prime(coef):
    """X(theta) = sum_k coef[k-1] cos(k theta), theta uniform: returns (k2, k3, k4). Exact (grid M > 4K)."""
    K = len(coef)
    M = 4 * K + 4
    th = 2 * np.pi * np.arange(M) / M
    X = np.zeros(M)
    for k, c in enumerate(coef, start=1):
        X += c * np.cos(k * th)
    m2 = np.mean(X ** 2); m3 = np.mean(X ** 3); m4 = np.mean(X ** 4)
    return m2, m3, m4 - 3 * m2 * m2


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, required=True)
    ap.add_argument('--smax', type=float, default=2.0)
    ap.add_argument('--sigmas', type=str, default='')
    ap.add_argument('--xmax', type=float, default=2.2e9, help='cap on q^s (sieve limit)')
    ap.add_argument('--seg', type=int, default=10_000_000)
    ap.add_argument('--out', type=str, required=True)
    ap.add_argument('--save-p', action='store_true')
    ap.add_argument('--ckpt-every', type=int, default=10)
    a = ap.parse_args()
    q = a.q
    assert is_prime(q) and q > 3
    logq = math.log(q)
    if a.sigmas:
        sig = np.array([float(x) for x in a.sigmas.split(',')])
    else:
        sig = np.round(np.arange(0.2, a.smax + 1e-9, 0.1), 10)
    sig = sig[q ** sig <= a.xmax * (1 + 1e-12)]
    ns = len(sig)
    X = int(math.floor(q ** sig.max()))          # all n < q^smax (fhat vanishes at u = s)
    g = primitive_root(q)
    L = dlog_table(q, g)
    root = math.isqrt(X)
    base = simple_sieve(root + 1)
    ck = a.out + '.ckpt.npz'
    t0 = time.time()
    if os.path.exists(ck):
        z = np.load(ck)
        assert int(z['q']) == q and np.allclose(z['sig'], sig) and int(z['X']) == X
        B = z['B']; D2 = z['D2']; D4 = z['D4']; nxt = int(z['nxt'])
        print(f'resume at {nxt}', flush=True)
    else:
        B = np.zeros((ns, q - 1))
        D2 = np.zeros(ns); D4 = np.zeros(ns)   # diagonal sums over primes p > sqrt(X) (only k = 1 matters there)
        nxt = 1
    # segments over [1, X]
    nseg_done = 0
    while nxt <= X:
        lo, hi = nxt, min(X + 1, nxt + a.seg)
        P = odd_primes_in(lo, hi, base)
        if lo <= 2 < hi:
            P = np.concatenate([np.array([2], dtype=np.int64), P])
        P = P[P != q]
        if len(P):
            lp = np.log(P.astype(np.float64))
            w0 = lp / np.sqrt(P.astype(np.float64))
            u = lp / logq
            l = L[P % q]
            big = P > root
            for j in range(ns):
                idx = np.searchsorted(u, sig[j])
                if idx == 0:
                    continue
                w = w0[:idx] * (1.0 - u[:idx] / sig[j])
                B[j] += np.bincount(l[:idx], weights=w, minlength=q - 1)
                bj = big[:idx]
                A = (2.0 / logq) * w[bj]
                D2[j] += np.sum(A * A) / 2.0
                D4[j] += -3.0 / 8.0 * np.sum(A ** 4)
        nxt = hi
        nseg_done += 1
        if nseg_done % a.ckpt_every == 0 and nxt <= X:
            np.savez(ck, q=q, sig=sig, X=X, B=B, D2=D2, D4=D4, nxt=nxt)
            print(f'  {nxt:.3e}/{X:.3e}  {time.time()-t0:.0f}s', flush=True)
    # prime powers p^k, k >= 2 (p <= sqrt(X)); and exact diagonal cumulants for the small primes
    Dsmall = np.zeros((ns, 3))
    for p in base:
        p = int(p)
        if p > root:
            break
        if p == q:
            continue
        lpp = math.log(p)
        K = int(math.floor(math.log(X) / lpp + 1e-12))
        while p ** (K + 1) <= X:
            K += 1
        while p ** K > X:
            K -= 1
        for j in range(ns):
            coef = []
            for k in range(1, K + 1):
                uk = k * lpp / logq
                fh = max(0.0, 1.0 - uk / sig[j])
                ak = lpp / p ** (k / 2.0) * fh
                coef.append(-(2.0 / logq) * ak)
                if k >= 2 and fh > 0:
                    B[j, L[pow(p, k, q)]] += ak
            while coef and coef[-1] == 0.0:
                coef.pop()
            if coef:
                k2, k3, k4 = diag_cumulants_one_prime(coef)
                Dsmall[j] += (k2, k3, k4)
    diag = np.zeros((ns, 3))
    diag[:, 0] = D2 + Dsmall[:, 0]
    diag[:, 1] = Dsmall[:, 1]
    diag[:, 2] = D4 + Dsmall[:, 2]
    # FFT per s
    k = np.arange(1, q - 1)
    even = (k % 2 == 0)
    psums = np.zeros((ns, 3, 9))
    Pall = np.zeros((ns, q - 2), dtype=np.float64) if a.save_p else None
    for j in range(ns):
        S = np.fft.ifft(B[j]) * (q - 1)        # S(k) = sum_l B[l] e^{+2 pi i k l/(q-1)}
        Pk = -(2.0 / logq) * S.real[1:q - 1]
        for c, m in enumerate((np.ones_like(even), even, ~even)):
            x = Pk[m]
            pw = np.ones_like(x)
            for e in range(9):
                psums[j, c, e] = pw.sum()
                pw = pw * x
        if Pall is not None:
            Pall[j] = Pk
    out = dict(q=q, g=g, sig=sig, X=X, psums=psums, diag=diag, logq=logq,
               runtime=time.time() - t0, fhat='tri')
    if Pall is not None:
        out['P'] = Pall
    np.savez(a.out, **out)
    if os.path.exists(ck):
        os.remove(ck)
    print(f'done q={q} g={g} X={X:.3e} ns={ns} time={time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
