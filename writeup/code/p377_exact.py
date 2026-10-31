# Second, separate check for writeup/377.md: exact counts by sieve (k <= 8). Needs ~300 MB of memory.
# Exact counts of solutions of variations #2 (n + 2*10^i prime for every digit position, all digits <= 7) and
# #3 (n - 2*10^i prime, all digits >= 2), by length k <= 8, base n composite or prime (counted separately).
import numpy as np, itertools
LIM = 3 * 10**8
s = np.ones(LIM, bool); s[:2] = False
for i in range(2, int(LIM**.5) + 1):
    if s[i]: s[i*i::i] = False
for var in (2, 3):
    for k in range(1, 9):
        digs = range(0, 8) if var == 2 else range(2, 10)
        lead = [d for d in digs if d >= 1]
        # all k-digit numbers with digits in digs (vectorized)
        n = np.array(lead, np.int64)
        for _ in range(k - 1): n = (n[:, None] * 10 + np.array(list(digs))[None, :]).ravel()
        ok = np.ones(len(n), bool)
        for i in range(k):
            m = n + (2 if var == 2 else -2) * 10**i
            ok &= s[m]
        sol = n[ok]
        print(f"var#{var} k={k}: {len(sol)} solutions ({int(s[sol].sum())} with prime base)", sol[-3:].tolist(), flush=True)
