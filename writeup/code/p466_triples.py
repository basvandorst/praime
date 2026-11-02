"""Puzzle 466: every m with phi(m-n) = phi(m) = phi(m+n), m <= N, for n = 1..NMAX. Usage: p466_triples.py 100000000 36 (15 s, 1.3 GB)."""
import numpy as np, sys
N = int(sys.argv[1])
phi = np.arange(N + 1, dtype=np.int32)
sv = np.ones(N + 1, dtype=np.bool_); sv[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
ps = np.nonzero(sv)[0]
del sv
for p in ps:
    p = int(p)
    v = phi[p::p]
    v -= v // p
for n in range(1, int(sys.argv[2]) + 1):
    a, b, c = phi[:N + 1 - 2 * n], phi[n:N + 1 - n], phi[2 * n:]
    ms = np.nonzero((a == b) & (b == c))[0] + n
    ms = ms[ms > n]
    print(n, len(ms), ms[:8].tolist(), flush=True)
