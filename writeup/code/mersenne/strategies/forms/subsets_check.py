"""Independent check of subsets.py: representations from PARI qfbsolve,
primes from sympy, plain Python.
"""
import json

import cypari2
from sympy import isprime, primerange

pari = cypari2.Pari()
N = 500_000


def solve(P, p):
    s = pari.qfbsolve(pari.Qfb(1, 0, P), p)
    if len(s) == 0:
        return None
    x, y = abs(int(s[0])), abs(int(s[1]))
    assert x * x + P * y * y == p
    return x, y


def main():
    ref = json.load(open("subsets.json"))
    primes = list(primerange(2, N + 1))
    q = [p for p in primes if p % 4 == 1]
    chk = {}
    e = {}
    for p in q:
        x, y = solve(1, p)
        e[p] = x if x % 2 == 0 else y
    r16 = {p: solve(16, p) is not None for p in q}
    r64 = {p: solve(64, p) is not None for p in q}
    chk["P=16 <=> 4|e"] = all(r16[p] == (e[p] % 4 == 0) for p in q)
    chk["P=16 <=> p=1 mod 8"] = all(r16[p] == (p % 8 == 1) for p in q)
    chk["P=64 <=> 8|e"] = all(r64[p] == (e[p] % 8 == 0) for p in q)
    chk["P=64 count"] = sum(r64.values()) == ref["P=64 count"]
    chk["P=16 count"] = sum(r16.values()) == ref["P=16 count"]
    q3 = [p for p in primes if p % 3 == 1]
    bq = {p: solve(3, p)[1] for p in q3}
    r12 = {p: solve(12, p) is not None for p in q3}
    r27 = {p: solve(27, p) is not None for p in q3}
    chk["P=12 <=> b even"] = all(r12[p] == (bq[p] % 2 == 0) for p in q3)
    chk["P=27 <=> 3|b"] = all(r27[p] == (bq[p] % 3 == 0) for p in q3)
    chk["P=27 count"] = sum(r27.values()) == ref["P=27 count"]
    run = 0
    while isprime(2 * run * run + 29):
        run += 1
    chk["2x^2+29 run"] = run == ref["2x^2+29 prime for x = 0..(run-1); run"]
    cls = sorted({p % 232 for p in primes if 232 % p and solve(58, p) is not None})
    chk["P=58 classes"] = cls == ref["P=58 classes mod 232"]
    # the other genus-mate: 2x^2 + 29y^2 represents the classes not in cls among (-58/p)=1
    other = sorted({p % 232 for p in primes if 232 % p and
                    len(pari.qfbsolve(pari.Qfb(2, 0, 29), p)) > 0})
    chk["2x^2+29y^2 classes disjoint from x^2+58y^2"] = not (set(other) & set(cls))
    chk["2x^2+29y^2 number of classes"] = len(other)
    g5 = [p for p in primes if 20 % p]
    chk["P=5 classes with (-5/p)=1"] = sorted({p % 20 for p in g5 if int(pari.kronecker(-5, p)) == 1}) \
        == ref["P=5 classes with (-5/p)=1"]
    chk["P=5 classes of x^2+5y^2"] = sorted({p % 20 for p in g5 if solve(5, p) is not None}) \
        == ref["P=5 classes of x^2+5y^2"]
    Q2 = pari.Qfb(2, 2, 3)
    chk["P=5 classes of 2x^2+2xy+3y^2"] = sorted({p % 20 for p in g5 if len(pari.qfbsolve(Q2, p)) > 0}) \
        == ref["P=5 classes of 2x^2+2xy+3y^2"]
    print(json.dumps(chk, indent=1))
    json.dump(chk, open("subsets_check.json", "w"), indent=1)


if __name__ == "__main__":
    main()
