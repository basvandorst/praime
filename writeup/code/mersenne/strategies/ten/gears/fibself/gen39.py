import json, re, os
HERE = os.path.dirname(os.path.abspath(__file__))
from sympy import primerange, isprime
R = json.load(open("ref39.json")); I = json.load(open("indices.json"))
L = [int(x) for x in open("fibpsp_1e+08.txt").read().split()]
pm1 = sum(1 for n in L if n % 5 in (1, 4)); pm2 = len(L) - pm1
g7 = sum(1 for n in L if n < 10**7 and n % 5 in (2, 3))
R5 = [int(x) for x in open(os.path.join(HERE, "..", "tour", "single_ring5_1e8.txt")).read().split("\n")[1:] if x.strip()]
S = set(L); r5ok = all(n in S and pow(2, n - 1, n) == 1 for n in R5)
pair = sum(1 for n in L if n % 5 in (2, 3) and pow(2, n - 1, n) == 1)
def alpha(m):
    a, b, k = 0, 1, 0
    while True:
        a, b, k = b, (a + b) % m, k + 1
        if a == 0: return k
P = list(primerange(2, 251)); A = {p: alpha(p) for p in P}
mers = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127]; tw = [p for p in P if isprime(p - 2) or isprime(p + 2)]
def best(T):
    bt, bc = 0, -1
    for t in range(1, 100001):
        c = sum(1 for p in T if t % A[p] == 0)
        if c > bc: bc, bt = c, t
    return bc, bt
b1, b2, b3 = best(P), best(mers), best(tw)
assert (b1, b2, b3) == ((26, 83160), (11, 27720), (19, 83160)) and r5ok and len(L) == 1460 and L[0] == 323 and g7 == 192 and pair == 0
f = lambda x: f"{x:,}"
recount = [
 ["Machine 6 (Fibonacci wheels): most primes at home 26 (tick 83,160), Mersenne exponents 11 (tick 27,720), twin primes 19 (tick 83,160)",
  f"the same: {b1[0]} at tick {f(b1[1])}, {b2[0]} at {f(b2[1])}, {b3[0]} at {f(b3[1])} (wheel p is home at tick t when the first Fibonacci number divisible by p, F(α(p)), has α(p) dividing t)"],
 ["Fibonacci wheel: first mistake 323, 1,460 mistakes below 10<sup>8</sup>",
  f"the same: {f(len(L))}, first {L[0]}; {f(pm1)} with n ≡ ±1 (mod 5) and {f(pm2)} with n ≡ ±2; the first 34 agree with OEIS A081264"],
 ["Every mistake of Pascal's wheel on 5 cells is a mistake of the Fibonacci and doubling wheels (90 below 10<sup>8</sup>)",
  f"the same, for all {len(R5)}"],
 ["Doubling and Fibonacci wheels together, n ≡ 2 or 3 (mod 5): no mistake below 2<sup>64</sup>",
  f"{pair} below 10<sup>8</sup>; the run to 2<sup>64</sup> uses the published list of base-2 pseudoprimes and was not repeated here"],
 ["Golden wheel (Pascal's 5 cells folded onto 2): 192 mistakes below 10<sup>7</sup>",
  f"the same, {g7}, but that run took only n ≡ 2 or 3 (mod 5), the numbers that do not fit 5 cells; the Fibonacci wheel makes exactly these {g7} mistakes there. The row on prime3 now says so"],
]
D = {"fibPrimeK": R["fib_prime_k"], "indices": I, "recount": recount}
src = open("prime39.src.html").read(); out = src.replace("__DATA__", json.dumps(D, separators=(",", ":"), ensure_ascii=False))
assert "__DATA__" not in out and "/*" not in out and "<!--" not in out
js = out.split("<script>")[1]; assert not re.search(r"(?<![:\\])//", re.sub(r"https?://[^\s\"'<>]+", "", js))
open("prime39.html", "w").write(out); print("prime39.html", len(out))
