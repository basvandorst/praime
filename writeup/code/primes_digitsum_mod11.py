import sys, time
import numpy as np
from math import sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve

isp = sieve(10**8)


def digit_sum(n):
    s = np.zeros_like(n); m = n.copy()
    while m.any():
        s += m % 10; m //= 10
    return s


print("E-F5e mechanism: 10 = -1 mod 11, so the alternating digit sum a(n) = n mod 11 and a(n) has the parity of the digit sum; P(even digit sum | n mod 11 = r) for all integers in [1e7, 1e8):")
n = np.arange(10**7, 10**8, 7)
s = digit_sum(n)
for r in range(11):
    m = n % 11 == r
    print(f"  r={r:2d}: P(even) = {(s[m] % 2 == 0).mean():.4f}", end="")
print()
print("primes avoid r = 0 (and multiples of 7, 13: 10^3 = -1 mod 7*11*13). Controls with more small primes removed:")
for lo, hi in ((10**5, 10**6), (10**6, 10**7), (10**7, 10**8)):
    pr = np.flatnonzero(isp[lo:hi]) + lo
    pe = (digit_sum(pr) % 2 == 0).mean()
    n = np.arange(lo, hi)
    line = f"  [{lo:.0e},{hi:.0e}): primes {pe:.5f} +- {sqrt(0.25/len(pr)):.5f}"
    for mod, label in ((30, "coprime to 30"), (2310, "coprime to 2310"), (30030, "coprime to 30030"), (9699690, "coprime to 19#")):
        keep = np.ones(len(n), dtype=bool)
        for q in (2, 3, 5, 7, 11, 13, 17, 19):
            if mod % q == 0:
                keep &= n % q != 0
        c = n[keep]
        ce = (digit_sum(c) % 2 == 0).mean()
        line += f" | {label} {ce:.5f} +- {sqrt(0.25/len(c)):.5f}"
    print(line, flush=True)
