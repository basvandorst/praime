import numpy as np, math
def count(a, b):
    """number of primes in [a, b), a even, plain segmented sieve"""
    s = np.ones(b - a, bool)
    r = int(math.isqrt(b)) + 1
    sp = np.ones(r + 1, bool); sp[:2] = False
    for i in range(2, int(r**0.5) + 1):
        if sp[i]: sp[i*i::i] = False
    for p in np.nonzero(sp)[0]:
        p = int(p); st = max(p*p, ((a + p - 1)//p)*p)
        s[st - a::p] = False
    return int(s.sum())
