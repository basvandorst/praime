"""Puzzle 1271: primes n with digits d_1 ... d_k such that (n - i) followed by the digit d_i is prime for i = 1..k
(91997 -> 919969, 919951, 919949, 919939, 919927).  Lists every solution below 10^7 (sieve to 10^8 for the
concatenations) and prints the heuristic expected number of k-digit solutions.
"""
import numpy as np, math
N = 10**7
s = np.ones(10 * N + 10, dtype=bool); s[:2] = False
for i in range(2, int(len(s)**0.5) + 1):
    if s[i]:
        s[i*i::i] = False
sol = {}
for p in np.flatnonzero(s[:N]).tolist():
    st = str(p)
    if any(c not in '1379' for c in st):       # the concatenations end in d_i, so every digit is 1, 3, 7 or 9
        continue
    k = len(st)
    if all(s[(p - i - 1) * 10 + int(st[i])] for i in range(k)):
        sol.setdefault(k, []).append(p)
assert sol == {1: [3, 7], 2: [19], 3: [773], 5: [91997]}, sol
print("solutions below 10^7 by number of digits:", sol)
# heuristic: 4^k strings over {1,3,7,9}; the string is prime with probability (5/2)/ln(10^k) (odd, prime to 5);
# each of the k concatenations is odd, prime to 5, of size 10^k: probability (5/2)/ln(10^k) each, independent
tot = 0
print("expected number of k-digit solutions:")
for k in range(2, 21):
    q = 2.5 / (k * math.log(10))
    e = 4**k * q * q**k
    tot += e if k >= 14 else 0
    print(f"  k={k:2d}: {e:.3g}" + ("   (found: %d)" % len(sol.get(k, [])) if k <= 7 else ""))
print(f"expected solutions with 14 or more digits (beyond the searched 1.8*10^13): {tot:.2g}")
