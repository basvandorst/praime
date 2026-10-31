"""Puzzle 938: soothsayer sequences. A k-digit N whose i-th digit equals omega(N+i) (I, distinct primes, A323083)
or Omega(N+i) (II, with multiplicity, A323084) for i = 1..k.
Checks every term on the page and the OEIS terms a(13)-a(17) of II (Resta 2019), then compares the count of terms
per digit length with the heuristic sum over k-digit N of prod_i P(digit_i = f(N+i)) (estimated by sampling)."""
import random
from sympy import factorint, primefactors
w = lambda n: len(primefactors(n)); W = lambda n: sum(factorint(n).values())
def ok(N, f): s = str(N); return all(int(c) == f(N + i) for i, c in enumerate(s, 1))
I = [0, 1, 12, 21, 22422, 24223, 33333, 34441524, 4242436235, 23443535352, 34462443242, 35256523324,
     4341535435353, 4645441523344, 5244526446515, 5335524234335]
II = [0, 1, 12, 21, 4224, 153426, 442451, 471614, 523291, 4336232, 474335342, 3624263478,
      36443455382, 244936365228, 452527642826, 593326437534, 4372566243537]
for N in I: assert ok(N, w), N
for N in II: assert ok(N, W), N
print("all", len(I), "terms of I and", len(II), "terms of II verified (II a(13)-a(17) from OEIS A323084)")
# heuristic: per position the chance that a digit matches f(n) is sum_d P(digit = d) P(f = d); digits ~ uniform 0..9
random.seed(938)
for name, f, seq in (("I", w, I), ("II", W, II)):
    row = []
    for k in range(5, 14):
        lo = 10 ** (k - 1); ns = [random.randrange(lo, 10 * lo) for _ in range(3000)]
        dist = [sum(f(n) == d for n in ns) / len(ns) for d in range(10)]
        p = sum(dist[1:]) / 10                       # leading digits 1..9 and others 0..9: about the same
        exp = 9 * lo * p ** k                        # 9*10^(k-1) candidates, k positions
        row.append(f"k={k}: expected {exp:.2f}, found {sum(len(str(N)) == k for N in seq)}")
    print(name, "; ".join(row))
