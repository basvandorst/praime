"""Independent reference for prime12.html (plain Python integers and fractions-free floats): after each cloud p
(primes in increasing order up to x) the number of dry numbers in [2, x] (no prime factor <= p) and the forecast
x prod_{q <= p} (1 - 1/q).   usage: ref12.py x out.json"""
import sys, json
from sympy import primerange
x = int(sys.argv[1])
wet = bytearray(x + 1); dry = x - 1; prod = 1.0; hist = []
for p in primerange(2, x + 1):
    for m in range(p, x + 1, p):
        if not wet[m]: wet[m] = 1; dry -= 1
    prod *= 1 - 1 / p
    hist.append([p, dry, x * prod])
assert hist[-1][1] == 0
json.dump(hist, open(sys.argv[2], 'w'))
for p, d, fc in hist:
    if p in (17, 31, 47, 67, 101, 251, 701, 997, 1999): print(p, d, round(fc, 1), round(d / fc, 4))
