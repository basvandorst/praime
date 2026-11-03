import sys
import numpy as np
from math import sqrt
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve

isp = sieve(2**27)


def ds2(n):
    s = np.zeros_like(n); m = n.copy()
    while m.any():
        s += m & 1; m >>= 1
    return s


print("E-F6b base 2 per bit length: P(even binary digit sum) for primes vs integers coprime to 6, 30, 2310, and the per-class table P(even | n mod 3 = r):")
for k in (20, 22, 24, 26):
    lo, hi = 2**k, 2**(k + 1)
    pr = np.flatnonzero(isp[lo:hi]) + lo
    n = np.arange(lo, hi)
    odd3 = (n % 2 == 1) & (n % 3 != 0)
    c6 = n[odd3]
    c2310 = n[odd3 & (n % 5 != 0) & (n % 7 != 0) & (n % 11 != 0)]
    c30030 = n[odd3 & (n % 5 != 0) & (n % 7 != 0) & (n % 11 != 0) & (n % 13 != 0)]
    pe = (ds2(pr) % 2 == 0).mean()
    print(f"  {k+1}-bit: primes {pe:.4f} +- {sqrt(0.25/len(pr)):.4f} | coprime to 6 {(ds2(c6)%2==0).mean():.4f} | to 2310 {(ds2(c2310)%2==0).mean():.4f} | to 30030 {(ds2(c30030)%2==0).mean():.4f}", flush=True)
n = np.arange(2**24, 2**25)
s = ds2(n) % 2 == 0
print("  25-bit integers: P(even | n mod 3 = r):", " ".join(f"r={r}: {s[n % 3 == r].mean():.4f}" for r in range(3)), "| P(even | n mod 5 = r):", " ".join(f"{s[n % 5 == r].mean():.4f}" for r in range(5)), "| mod 7:", " ".join(f"{s[n % 7 == r].mean():.4f}" for r in range(7)))
