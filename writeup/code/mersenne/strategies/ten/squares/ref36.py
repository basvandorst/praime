"""Independent reference for prime36.html (prime squares), plain integers and sympy.

Grid of N x N dots (x, y = 0..N-1), point zero top left. Numbering:
  rows:   n = x + N*y
  shells: k = max(x, y); n = k^2 + x if y == k, else k^2 + k + (k - y)   (shell k holds k^2 .. k^2 + 2k)
For every prime p < N the square from point zero to (p, p); its outline is the set of dots with
max(x, y) = p. A dot on an outline is 'crossed'.
Prints, per prime p, the primes on the outline (top edge y = p, right edge x = p), and totals.
"""
import json
import sys
from sympy import isprime, primerange


def number(x, y, N, mode):
    if mode == "rows":
        return x + N * y
    k = max(x, y)
    return k * k + x if y == k else k * k + k + (k - y)


def sizes(N, which):
    if which == "primes":
        return list(primerange(2, N))
    if which == "fakes":
        return [s for s in range(3, N) if s % 2 and s % 5 and not isprime(s)]
    return list(range(2, N, 2))


def run(N, mode, which="primes"):
    ps = sizes(N, which)
    per = {}
    crossed = set()
    for p in ps:
        top = [number(x, p, N, mode) for x in range(p + 1)]
        right = [number(p, y, N, mode) for y in range(p)]
        per[p] = [sum(map(isprime, top)), sum(map(isprime, right))]
        crossed.update((x, p) for x in range(p + 1))
        crossed.update((p, y) for y in range(p))
    allp = sum(isprime(number(x, y, N, mode)) for x in range(N) for y in range(N))
    cp = sum(isprime(number(x, y, N, mode)) for x, y in crossed)
    return {"per": per, "dots": N * N, "primes": allp, "crossed": len(crossed), "crossed_primes": cp}


if __name__ == "__main__":
    out = {}
    for N in (20, 50, 100):
        for mode, which in (("rows", "primes"), ("shells", "primes"), ("rows", "fakes"), ("shells", "fakes"), ("rows", "evens"), ("shells", "evens")):
            r = run(N, mode, which)
            out["%s_%s_%d" % (mode, which, N)] = r
            print("%s, squares %s, N=%d: %d dots, %d prime dots (%.4f); %d crossed dots, %d of them prime (%.4f); off the outlines %.4f" % (
                mode, which, N, r["dots"], r["primes"], r["primes"] / r["dots"], r["crossed"], r["crossed_primes"],
                r["crossed_primes"] / r["crossed"], (r["primes"] - r["crossed_primes"]) / (r["dots"] - r["crossed"])))
    r = out["rows_primes_100"]["per"]
    print("rows N=100, right edge x = p is p, p+100, ...: primes per p:", {p: v[1] for p, v in r.items()})
    print("rows N=100, top edge y = p is 100p .. 101p: primes per p:", {p: v[0] for p, v in r.items()})
    s = out["shells_primes_100"]["per"]
    print("shells N=100, outline p holds p^2 .. p^2+2p: primes per p:", {p: v[0] + v[1] for p, v in s.items()})
    json.dump(out, open("ref36.json", "w"))
