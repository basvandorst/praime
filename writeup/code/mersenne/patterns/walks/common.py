"""Shared sieve and series for the walk / turtle / row layouts of prime2.html.

All series are indexed by n (index 0 unused or 0), exactly as the page computes them:
  M(n)    = sum_{m<=n} mu(m)
  L(n)    = sum_{m<=n} lambda(m),  lambda = (-1)^Omega
  R(n)    = #{p<=n, p=3 mod 4} - #{p<=n, p=1 mod 4}
  turtle  = start at (0,0) for n=1; the step into n (n>=2) has heading pi(n-1) mod 4,
            DX=[1,0,-1,0], DY=[0,-1,0,1] (canvas y points down).
"""
import numpy as np

VIS = "/tmp/claude-1000/-home-bas-praime/d6f1eef8-404d-4ff2-8eca-89d091b2ed2f/scratchpad/mexp/vis"
N = 500000
EXPONENTS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
             4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503,
             132049, 216091]


def sieve(n):
    """smallest prime factor spf[0..n] (spf[1]=1), Omega, mu, prime mask."""
    spf = np.zeros(n + 1, dtype=np.int32)
    r = int(n ** 0.5)
    for p in range(2, r + 1):
        if spf[p] == 0:
            idx = np.arange(p * p, n + 1, p)
            free = spf[idx] == 0
            spf[idx[free]] = p
    ar = np.arange(n + 1, dtype=np.int32)
    z = spf == 0
    spf[z] = ar[z]
    spf[0] = 0
    spf[1] = 1
    isp = np.zeros(n + 1, dtype=bool)
    isp[2:] = spf[2:] == ar[2:]
    primes = np.nonzero(isp)[0]
    # Omega: add 1 for every prime power p^k dividing n
    Om = np.zeros(n + 1, dtype=np.int16)
    sqf = np.ones(n + 1, dtype=bool)
    for p in primes:
        q = int(p)
        k = 1
        while q <= n:
            Om[q::q] += 1
            if k == 2:
                sqf[q::q] = False
            q *= int(p)
            k += 1
    lam = np.where(Om % 2 == 0, 1, -1).astype(np.int8)
    lam[0] = 0
    mu = np.where(sqf, lam, 0).astype(np.int8)
    mu[0] = 0
    return spf, Om, mu, lam, isp, primes


def series(n=N):
    spf, Om, mu, lam, isp, primes = sieve(n)
    M = np.cumsum(mu.astype(np.int64))
    L = np.cumsum(lam.astype(np.int64))
    ar = np.arange(n + 1)
    step = np.zeros(n + 1, dtype=np.int64)
    step[isp & (ar % 4 == 3)] = 1
    step[isp & (ar % 4 == 1)] = -1
    R = np.cumsum(step)
    pi = np.cumsum(isp.astype(np.int64))
    # turtle: heading of the step into n is pi(n-1) mod 4
    head = np.zeros(n + 1, dtype=np.int64)
    head[2:] = pi[1:n] % 4
    DX = np.array([1, 0, -1, 0]); DY = np.array([0, -1, 0, 1])
    dx = np.zeros(n + 1, dtype=np.int64); dy = np.zeros(n + 1, dtype=np.int64)
    dx[2:] = DX[head[2:]]; dy[2:] = DY[head[2:]]
    X = np.cumsum(dx); Y = np.cumsum(dy)
    return dict(spf=spf, Om=Om, mu=mu, lam=lam, isp=isp, primes=primes, M=M, L=L, R=R,
                pi=pi, head=head, X=X, Y=Y)


def load_pos(layout):
    fx = np.fromfile(f"{VIS}/data/{layout}.fx.bin", dtype="<f4")
    fy = np.fromfile(f"{VIS}/data/{layout}.fy.bin", dtype="<f4")
    return fx, fy
