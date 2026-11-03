import sys, time
import numpy as np
from math import log, sqrt, pi, atan2, gcd
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, primes_upto, hl_constant, lam_array, C2, GAMMA

t0 = time.time()
X = 10**8
isp = sieve(2 * X + 10)
pr = np.flatnonzero(isp[:X + 1])
pr = np.flatnonzero(isp[:X + 1])
print("E-C17 Hardy-Littlewood's second conjecture pi(x+y) <= pi(x) + pi(y) for y <= 3000 and x <= 1e8 (false far out, Hensley-Richards):")
pic = np.cumsum(isp[:X + 1])
worst = (0, 0, 0)
for y in (100, 300, 1000, 3000):
    mx = (pic[y + 2:X + 1] - pic[2:X + 1 - y]).max()
    print(f"  y={y}: max over x >= 2 of pi(x+y) - pi(x) = {mx}, pi(y) = {pic[y]}, margin {pic[y] - mx}")
print(f"[{time.time()-t0:.0f} s]")
