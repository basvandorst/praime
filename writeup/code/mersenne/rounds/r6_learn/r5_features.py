"""Frozen feature definitions (round-3 near-misses) and the broad feature family.

All functions take a numpy int64 array of odd primes n and return a float array
(nan = outside the feature's domain). Class tests are in CLASSES.
"""
import math
import numpy as np


def mod7(n):
    return (n % 7 == 1).astype(float)


def angle_x2_2y2(n):
    """For n = 1, 3 (mod 8): n = x^2 + 2y^2 (x, y > 0, unique for prime n);
    angle of (x, y*sqrt2) divided by pi/2, in [0, 1]. nan otherwise."""
    import sympy
    out = np.full(len(n), np.nan)
    for i, p in enumerate(n.tolist()):
        if p % 8 not in (1, 3):
            continue
        # Cornacchia: r^2 = -2 mod p, Euclid on (p, r) down to the first remainder < sqrt p
        r0 = sympy.sqrt_mod(-2 % p, p)
        if r0 is not None:
            a_, b_ = p, r0
            while b_ * b_ >= p:
                a_, b_ = b_, a_ % b_
            x = b_
            y2 = (p - x * x) // 2
            y = math.isqrt(y2) if y2 >= 0 else -1
            if y > 0 and x * x + 2 * y * y == p:
                out[i] = math.atan2(math.sqrt(2) * y, x) / (math.pi / 2)
                continue
        y = 1
        while 2 * y * y < p:
            x2 = p - 2 * y * y
            x = math.isqrt(x2)
            if x * x == x2:
                out[i] = math.atan2(math.sqrt(2) * y, x) / (math.pi / 2)
                break
            y += 1
    return out


def frac_inv2pi(n):
    x = n.astype(np.float64) / (2 * np.pi)       # n <= 1.4e8: float64 error ~1e-8, fine
    return x - np.floor(x)


def cos_n(n):
    return np.cos(n.astype(np.float64))


def collatz_glide(n):
    """Standard map (n/2 or 3n+1); first step t with value < n."""
    out = np.zeros(len(n))
    for i, k in enumerate(n.tolist()):
        v, t = k, 0
        while True:
            v = v // 2 if v % 2 == 0 else 3 * v + 1
            t += 1
            if v < k:
                break
        out[i] = t
    return out


def cf_maxq_norm(n):
    """Largest partial quotient of sqrt(n) inside the period (the final 2*a0 excluded;
    0 for period length 1), divided by sqrt(n)."""
    out = np.zeros(len(n))
    for i, k in enumerate(n.tolist()):
        a0 = math.isqrt(k)
        m, d, a = 0, 1, a0
        mx = 0
        while True:
            m = d * a - m
            d = (k - m * m) // d
            a = (a0 + m) // d
            if a == 2 * a0:
                break
            mx = max(mx, a)
        out[i] = mx / math.sqrt(k)
    return out


# (name, function, class test on the feature value, domain note)
FROZEN = [
    ("H1 n=1 mod 7", mod7, lambda v: v == 1),
    ("H2 x^2+2y^2 angle >0.75 (n=1,3 mod 8)", angle_x2_2y2, lambda v: v > 0.75),
    ("H3 frac(n/2pi)<=0.103", frac_inv2pi, lambda v: v <= 0.103),
    ("H4 cos n>0.803", cos_n, lambda v: v > 0.803),
    ("H5 Collatz glide>21", collatz_glide, lambda v: v > 21),
    ("H6 cf sqrt maxq/sqrt n<=0.397", cf_maxq_norm, lambda v: v <= 0.397),
]


def digit_features(n):
    """Simple digit statistics; returns dict name -> boolean array (all defined)."""
    s = [str(int(x)) for x in n]
    bstr = [bin(int(x))[2:] for x in n]
    dsum = np.array([sum(map(int, t)) for t in s], float)
    pop = np.array([t.count("1") for t in bstr], float)
    runs = np.array([1 + sum(t[i] != t[i + 1] for i in range(len(t) - 1)) for t in bstr], float)
    ndist = np.array([len(set(t)) for t in s], float)
    lead = np.array([int(t[0]) for t in s])
    blen = np.array([len(t) for t in bstr])
    # size-free versions: z-scores against the binomial/uniform expectation for that length
    dlen = np.array([len(t) for t in s], float)
    dsum_z = (dsum - 4.5 * dlen) / np.sqrt(8.25 * dlen)
    pop_z = (pop - 1 - (blen - 1) / 2) / np.sqrt(np.maximum(blen - 1, 1) / 4)
    runs_z = (runs - 1 - (blen - 1) / 2) / np.sqrt(np.maximum(blen - 1, 1) / 4)
    rev10 = np.array([int(t[::-1]) for t in s])
    rev2 = np.array([int(t[::-1], 2) for t in bstr])
    import sympy
    f = {}
    f["dsum_z<=-1"] = dsum_z <= -1
    f["dsum_z>=1"] = dsum_z >= 1
    f["popcount_z<=-1"] = pop_z <= -1
    f["popcount_z>=1"] = pop_z >= 1
    f["popcount odd"] = pop % 2 == 1
    f["binary runs_z<=-1"] = runs_z <= -1
    f["binary runs_z>=1"] = runs_z >= 1
    f["distinct dec digits = len"] = ndist == dlen
    f["has dec digit 0"] = np.array(["0" in t for t in s])
    f["all dec digits odd"] = np.array([all(c in "13579" for c in t) for t in s])
    f["dec palindrome"] = np.array([t == t[::-1] for t in s])
    f["bin palindrome"] = np.array([t == t[::-1] for t in bstr])
    f["dec reversal prime"] = np.array([sympy.isprime(int(x)) for x in rev10])
    f["bin reversal prime"] = np.array([sympy.isprime(int(x)) for x in rev2])
    f["binary length odd"] = blen % 2 == 1
    f["dec length odd"] = dlen % 2 == 1
    for d in range(1, 10):
        f[f"leading dec digit {d}"] = lead == d
    f["leading two bits 10"] = np.array([t[:2] == "10" for t in bstr])
    f["trailing bits 1111"] = (n % 16) == 15
    f["n+2 prime"] = np.array([sympy.isprime(int(x) + 2) for x in n])
    f["n-2 prime"] = np.array([sympy.isprime(int(x) - 2) for x in n])
    f["2n+1 prime"] = np.array([sympy.isprime(2 * int(x) + 1) for x in n])
    f["2n-1 prime"] = np.array([sympy.isprime(2 * int(x) - 1) for x in n])
    return f
