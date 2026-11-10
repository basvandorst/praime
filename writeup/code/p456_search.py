# Puzzle 456: all ways to write M = 2^e - 1 as a sum of k >= 2 consecutive primes, for the Mersenne primes with e <= 31.
# Prefix sums C_j = p_1 + ... + p_j; a window is a pair i < j with C_j - C_i = M.  Two segmented prime streams
# (lead and lag) are run side by side; memory stays below 300 MB.
# usage: (ulimit -v 1500000; p456_search.py [e ...])      default: 2 3 5 7 13 17 19 31
import sys, numpy as np
SEG = 1 << 24
def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)
def stream(limit):
    """yields (primes, prefix sums) for the primes <= limit, segment by segment"""
    base = small_primes(int(limit ** .5) + 2); total = 0
    for lo in range(0, limit + 1, SEG):
        hi = min(lo + SEG, limit + 1); w = np.ones(hi - lo, dtype=bool)
        if lo == 0: w[:2] = False
        for p in base:
            p = int(p)
            if p * p >= hi: break
            w[max(p * p, (lo + p - 1) // p * p) - lo::p] = False
        P = np.flatnonzero(w).astype(np.int64) + lo
        if len(P) == 0: continue
        C = np.cumsum(P) + total; total = int(C[-1])
        yield P, C
def windows(M):
    out = []
    lag = stream(M); LP = np.array([0], dtype=np.int64); LC = np.array([0], dtype=np.int64)   # C_0 = 0 (window from 2)
    for P, C in stream(min(M, M // 2 + 2000)):                  # k >= 2: the largest prime is below M/2 + (a gap); checked: gaps < 2000 here
        T = C - M
        T = T[T >= 0]
        if len(T) == 0: continue
        while LC[-1] < T[-1]:
            lp, lc = next(lag); LP = np.concatenate([LP, lp]); LC = np.concatenate([LC, lc])
        keep = np.searchsorted(LC, T[0]); LP, LC = LP[keep:], LC[keep:]
        pos = np.searchsorted(LC, T); pos[pos >= len(LC)] = len(LC) - 1
        for t in np.flatnonzero(LC[pos] == T):
            j_prime = int(P[len(P) - len(T) + t]); last_before = int(LP[pos[t]])   # window = primes in (last_before, j_prime]
            if j_prime < M: out.append((last_before, j_prime))            # j_prime = M would be the single prime M
    return out
if __name__ == "__main__":
    for e in [int(x) for x in sys.argv[1:]] or [2, 3, 5, 7, 13, 17, 19, 31]:
        M = 2 ** e - 1
        print(f"2^{e} - 1 = {M}: {len(windows(M))} windows (largest prime before the window, last prime):", windows(M), flush=True)
