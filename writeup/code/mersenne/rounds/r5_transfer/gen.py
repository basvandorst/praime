#!/usr/bin/env python3
"""Generate PRP data for families with the Mersenne structure.

Family ('m', b): R = (b^n - 1)/(b - 1);  family ('p', b): R = (b^n + 1)/(b + 1)
(('p', 2) is the Wagstaff family). b is never a perfect power. Only odd prime n.
For every odd prime n <= N_b = floor(BITS / log2 b):
  * n | b -/+ 1  -> n divides R: composite, status 'ndiv'
  * sieve: q = 2kn + 1, k = 1..K (every prime factor of R is = 1 mod 2n), skipping q with a
    factor <= 23; q | R  <=>  b^n = +/-1 (mod q) and gcd(q, b -/+ 1) = 1 and q < R
    -> status 'factor' with the smallest such q (q itself may be composite; then its prime
    factors, also = 1 mod 2n, divide R)
  * otherwise a Fermat test base a (a = 3, or 2 when b = 3) -> 'fermat' if it fails,
    else gmpy2.is_prime (BPSW + MR) -> 'PRP' or 'mr'.
Output: one JSON line per family in data.jsonl (append; families already present are skipped).
"""
import json, math, os, sys, time
import numpy as np
import gmpy2
import sympy

BITS = 10000
K = 100000
OUT = os.environ.get("GEN_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "data.jsonl"))
SMALL = [3, 5, 7, 11, 13, 17, 19, 23]


def is_perfect_power(b):
    return any(round(b ** (1.0 / e)) ** e == b for e in range(2, 8))


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def vpowmod(base, e, q):
    """base^e mod q, elementwise; q < 2^31 so products fit in uint64."""
    r = np.ones_like(q)
    x = base % q
    while e:
        if e & 1:
            r = (r * x) % q
        x = (x * x) % q
        e >>= 1
    return r


def sieve_n(kind, b, n, Rexact):
    k = np.arange(1, K + 1, dtype=np.uint64)
    q = np.uint64(2 * n) * k + np.uint64(1)
    keep = np.ones(len(q), bool)
    for p in SMALL:
        keep &= (q % np.uint64(p) != 0) | (q == np.uint64(p))
    q = q[keep]
    r = vpowmod(np.full(len(q), b, dtype=np.uint64), n, q)
    target = np.ones_like(q) if kind == "m" else q - np.uint64(1)
    hit = r == target
    if not hit.any():
        return None
    g = b - 1 if kind == "m" else b + 1
    for qq in q[hit]:
        qq = int(qq)
        if math.gcd(qq, g) == 1 and (Rexact is None or qq < Rexact):
            return qq
    return None


def run_family(kind, b):
    N = int(BITS / math.log2(b))
    g = b - 1 if kind == "m" else b + 1
    a = 2 if b % 3 == 0 and is_perfect_power_of(b, 3) else 3
    recs = []
    t0 = time.time()
    nf = 0
    for n in primes_upto(N):
        n = int(n)
        if n == 2:
            continue
        B = gmpy2.mpz(b)
        R = (B ** n - 1) // (b - 1) if kind == "m" else (B ** n + 1) // (b + 1)
        if g % n == 0 and R != n:
            recs.append((n, "ndiv", n)); continue
        Rex = int(R) if n * math.log2(b) < 80 else None
        q = sieve_n(kind, b, n, Rex)
        if q is not None:
            recs.append((n, "factor", q)); continue
        nf += 1
        if R < 2 ** 64:                       # proof range: decide exactly
            recs.append((n, "PRP" if sympy.isprime(int(R)) else "fermat", a)); continue
        if gmpy2.powmod(a, R - 1, R) != 1:
            recs.append((n, "fermat", a)); continue
        recs.append((n, "PRP" if gmpy2.is_prime(R) else "mr", None))
    return {"kind": kind, "b": b, "N": N, "K": K, "bits": BITS, "fermat_base": a,
            "n_fermat_tests": nf, "secs": round(time.time() - t0, 1), "recs": recs}


def is_perfect_power_of(b, p):
    while b % p == 0:
        b //= p
    return b == 1


def families():
    fam = []
    for b in range(3, 51):
        if not is_perfect_power(b):
            fam.append(("m", b))
    for b in range(2, 51):
        if not is_perfect_power(b):
            fam.append(("p", b))
    return fam


if __name__ == "__main__":
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            d = json.loads(line); done.add((d["kind"], d["b"]))
    fams = families()
    if len(sys.argv) > 1:
        fams = [f for f in fams if f"{f[0]}{f[1]}" in sys.argv[1:]]
    for kind, b in fams:
        if (kind, b) in done:
            continue
        d = run_family(kind, b)
        prp = [r[0] for r in d["recs"] if r[1] == "PRP"]
        print(f"{kind}{b}: N={d['N']} n={len(d['recs'])} fermat={d['n_fermat_tests']} "
              f"PRP={len(prp)} {prp} ({d['secs']}s)", flush=True)
        with open(OUT, "a") as f:
            f.write(json.dumps(d) + "\n")
