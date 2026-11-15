#!/usr/bin/env python3
"""Independent recheck of data.jsonl with plain Python integers (no gmpy2, no numpy):
  * every 'factor' record: R computed exactly, R % q == 0 and 1 < q < R
  * every 'ndiv' record: R % n == 0
  * every PRP with log2 R <= 4000, plus 40 random larger ones: sympy.isprime (a proof) when
    R < 2^64, else Fermat bases 5 (7 for b = 5) and 11, other than the generator's base
  * 150 random 'fermat' records: not prime (sympy below 2^64, else a failed Fermat test)
  * every record with R < 2^64: status agrees with sympy.isprime
  * completeness: the records of each family are exactly the odd primes <= N_b (sympy.primerange)
"""
import json, math, os, random
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
random.seed(7)
fams = [json.loads(l) for l in open(os.path.join(HERE, "data.jsonl"))]


def base(b):
    return 7 if b == 5 else 5          # never a power of b (b = 5 is the only power of 5 used)


def R(kind, b, n):
    return (b ** n - 1) // (b - 1) if kind == "m" else (b ** n + 1) // (b + 1)


nf = nd = nprp = nfer = 0; bad = []
big_prp, fer = [], []
for d in fams:
    kind, b, N = d["kind"], d["b"], d["N"]
    ns = [r[0] for r in d["recs"]]
    if ns != list(sympy.primerange(3, N + 1)):
        bad.append(("incomplete", kind, b))
    for n, st, x in d["recs"]:
        if st == "factor":
            v = R(kind, b, n); nf += 1
            if not (1 < x < v and v % x == 0):
                bad.append(("factor", kind, b, n, x))
        elif st == "ndiv":
            nd += 1
            if R(kind, b, n) % n:
                bad.append(("ndiv", kind, b, n))
        elif st == "PRP":
            if n * math.log2(b) <= 4000:
                nprp += 1
                v = R(kind, b, n)
                if v < 2 ** 64:
                    if not sympy.isprime(v):
                        bad.append(("prp<2^64", kind, b, n))
                elif not (pow(base(b), v - 1, v) == 1 and pow(11, v - 1, v) == 1):
                    bad.append(("prp", kind, b, n))
            else:
                big_prp.append((kind, b, n))
        elif st == "fermat":
            fer.append((kind, b, n))
for kind, b, n in random.sample(big_prp, min(40, len(big_prp))):
    v = R(kind, b, n); nprp += 1
    if not (pow(base(b), v - 1, v) == 1 and pow(11, v - 1, v) == 1):
        bad.append(("prp", kind, b, n))
for kind, b, n in random.sample(fer, min(150, len(fer))):
    v = R(kind, b, n); nfer += 1
    if (sympy.isprime(v) if v < 2 ** 64 else pow(base(b), v - 1, v) == 1):
        bad.append(("fermat", kind, b, n))
nsmall = 0
for d in fams:
    for n, st, x in d["recs"]:
        v = R(d["kind"], d["b"], n)
        if v < 2 ** 64:
            nsmall += 1
            if sympy.isprime(v) != (st == "PRP"):
                bad.append(("small", d["kind"], d["b"], n, st))
print(f"records with R < 2^64 proved with sympy: {nsmall}")
print(f"families {len(fams)}; factors checked {nf}; ndiv {nd}; PRPs rechecked {nprp} "
      f"(of {nprp - min(40, len(big_prp)) + len(big_prp)}); fermat composites rechecked {nfer}")
print("problems:", bad if bad else "none")
