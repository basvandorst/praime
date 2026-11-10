# Puzzles 655 / 498, checks.
# (1) Completeness of the Pell machinery: p655_pell.py --raw (no sieve, no primality) must list exactly the integer
#     solutions of 2p^2 + 2tp + 4C = square with 3 <= p < 10^7, half-gaps a < b < d <= 12, t = a + b + d odd, found
#     here by plain enumeration of p.
# (2) Every solution in SOL is recomputed: p, q, r, s proved prime (PARI isprime), nothing prime in between,
#     the form a perfect square (math.isqrt), and the Pell identity of the write-up.
# Run from the repo root: .venv/bin/python3 writeup/code/p655_verify.py
import subprocess, sys, math
import numpy as np
out = subprocess.run([sys.executable, 'writeup/code/p655_pell.py', '7', '12', '3', '--raw'], capture_output=True, text=True).stdout
pell = {eval(l) for l in out.splitlines() if l.startswith('    (')}
brute = set(); p = np.arange(3, 10**7, dtype=np.int64)
for d in range(3, 13):
    for b in range(2, d):
        for a in range(1, b):
            t = a + b + d
            if t % 2 == 0: continue
            for typ, C in ((1, b * d), (2, a * b), (3, a * d)):
                N = 2 * p * p + 2 * t * p + 4 * C
                r = np.rint(np.sqrt(N.astype(np.float64))).astype(np.int64)
                for i in np.flatnonzero(r * r == N): brute.add((int(p[i]), typ, a, b, d))
assert pell == brute, (len(pell), len(brute), sorted(pell ^ brute)[:5])
print(f"(1) Pell orbits = plain enumeration: {len(brute)} integer solutions with p < 10^7, d <= 12")
import cypari2
pari = cypari2.Pari()
SOL = [(203233, 1),
       (1150284232069964804556646431054740232316097446258616188423, 1),
       (11676583373868016849349945425661876543184231432061070831984678948459667093277, 2)]
for p, typ in SOL:
    # the next three primes after p: every number in between fails a probable-prime test (so it is composite),
    # and p, q, r, s are proved prime by PARI's isprime (APR-CL / ECPP, a proof at any size)
    primes = [p]; x = p
    while len(primes) < 4:
        x += 1
        if pari.ispseudoprime(x): primes.append(x)
    assert all(pari.isprime(x) for x in primes)
    p, q, r, s = primes
    N = {1: p * q + r * s, 2: p * s + q * r, 3: p * r + q * s}[typ]
    assert math.isqrt(N) ** 2 == N
    print(f"(2) type {typ}: p = {p} ({len(str(p))} digits, proved prime), q = p + {q - p}, r = p + {r - p}, s = p + {s - p} consecutive proved primes: {['', 'pq + rs', 'ps + qr', 'pr + qs'][typ]} = {math.isqrt(N)}^2")
    a, b, d = (q - p) // 2, (r - p) // 2, (s - p) // 2; t = a + b + d; C = {1: b * d, 2: a * b, 3: a * d}[typ]
    U = 2 * p + t; z = math.isqrt(N) // 2; k = t * t - 8 * C
    assert t % 2 == 1 and U * U - 8 * z * z == k
    print(f"    half-gaps {a}, {b}, {d}: t = {t}, k = t^2 - 8C = {k}, U = 2p + t, U^2 - 8 z^2 = k holds; ln p = {math.log(p):.1f}, (s - p)/ln p = {(s - p) / math.log(p):.2f}")
print("all checks passed")
