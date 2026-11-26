"""Natural-weight data: for N prime, weight k in {2,4}, every prime p <= YMAX*N^2 (p != N):
   s(p) = sum over newforms f in S_k(Gamma0(N)) of lambda_f(p)^2,  lambda_f(p) = a_f(p)/p^{(k-1)/2},
so that E_nat[lambda_f(p^2)] = s(p)/dim - 1.  Method: trace.Tracer (Eichler-Selberg, validated against PARI) and the
normalized Gram identity  s(p) = w^T G^{-1} w,  G_{mm'} = sum_f lambda_f(m) lambda_f(m'),  w_m = sum_f lambda_f(p) lambda_f(m),
m in S = {1..dim} (rank checked).  Usage: gen_data.py CODEDIR OUTDIR K YMAX BUDGET N1 N2 ...  (skips levels already saved)."""
import sys, os, time
sys.path.insert(0, sys.argv[1])
import numpy as np
from math import gcd
from sympy import primerange, divisors
from hurwitz import hurwitz12
from trace import Tracer

outdir, k, ymax, budget = sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
levels = [int(x) for x in sys.argv[6:]]
todo = [N for N in levels if not os.path.exists(os.path.join(outdir, "data_k%d_N%d.npz" % (k, N)))]
if not todo:
    print("nothing to do"); sys.exit(0)
t0 = time.time()
# H table large enough for the largest level of this run: n <= B*dim, dim <= (k-1)(N+1)/12 + 2
Nmax = max(todo)
B = int(ymax * Nmax * Nmax)
dmax = (k - 1) * (Nmax + 1) // 12 + 2
H = hurwitz12(4 * B * dmax + 16)
print("H table to %d in %.1fs" % (len(H) - 1, time.time() - t0), flush=True)
sig = {}
for N in todo:
    if time.time() - t0 > budget:
        break
    t1 = time.time()
    T = Tracer(N, k, H)
    S = T.basis()
    dim = T.dim
    Sa = np.array(S, dtype=float)
    nrm = Sa ** ((k - 1) / 2.0)
    G = np.array(T.Gint, dtype=float) / np.outer(nrm, nrm)
    Ginv = np.linalg.inv(G)
    sig_m = {m: sum(d ** (k - 1) + (m // d) ** (k - 1) for d in divisors(m)) for m in S}
    P = [p for p in primerange(2, int(ymax * N * N) + 1) if p != N]
    out = np.empty(len(P))
    for i, p in enumerate(P):
        w = np.empty(dim)
        for j, m in enumerate(S):
            if p > m:
                n = p * m
                # fast trace for n = p*m, p prime > m: not a square, A3 = -2 sigma_{k-1}(m), A4 = (1+p) sigma(m) (k=2)
                if n in T.cache:
                    v = T.cache[n]
                else:
                    from math import isqrt
                    r = isqrt(4 * n); tmax = r if r * r < 4 * n else r - 1
                    t = np.arange(-tmax, tmax + 1, dtype=np.int64)
                    D = 4 * n - t * t
                    term = H[D].astype(np.int64) * (1 + T.leg[(-D) % N])
                    msk = (D % (N * N) == 0)
                    if msk.any():
                        term[msk] += N * H[D[msk] // (N * N)]
                    S2 = int(term.sum()) if k == 2 else int(((t * t - n) * term).sum())
                    A3 = -sig_m[m]
                    A4 = (1 + p) * sum(divisors(m)) if k == 2 else 0
                    tot = -S2 + 24 * (A3 + A4)
                    assert tot % 24 == 0
                    v = tot // 24
                w[j] = v / (p * m) ** ((k - 1) / 2.0)
            else:
                w[j] = T.pair(p, m) / (p * m) ** ((k - 1) / 2.0)
        out[i] = w @ Ginv @ w
    np.savez(os.path.join(outdir, "data_k%d_N%d.npz" % (k, N)), p=np.array(P, dtype=np.int64), s=out, dim=dim, N=N, k=k,
             S=np.array(S), cond=np.linalg.cond(G))
    print("N=%d k=%d dim=%d primes=%d cond(G)=%.1e  %.1fs" % (N, k, dim, len(P), np.linalg.cond(G), time.time() - t1), flush=True)
print("total %.1fs" % (time.time() - t0))
