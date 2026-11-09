# Puzzle 1049: f(p) = sum of d^d over the digits of p; chain = successive zero-free primes.
# (1) page values n = 1..3 and n = 4; (2) minimality of n = 4 recomputed with sympy over all digit multisets of
# <= 12 digits; (3) n = 5: the 321-digit p of data/p1049_n5_strict.txt: PARI isprime proof, its chain, the digit
# count D(q) of every candidate q recomputed by an independent exact coin-change method, and the 318 smaller arrangements composite.
# Run from the repo root: .venv/bin/python3 writeup/code/p1049_verify.py
import itertools, re
import numpy as np
from sympy import isprime
import cypari2
pari = cypari2.Pari(); pari.allocatemem(600 * 2**20)
POW = [d**d for d in range(10)]
def f(x): return sum(POW[int(c)] for c in str(x))
def zf_prime(x): return '0' not in str(x) and isprime(x)
def chain(p):
    k = 0
    while True:
        q = f(p)
        if not zf_prime(q): return k
        k += 1; p = q

for p, n in ((11, 1), (1433, 2), (122527, 3), (112222268681, 4)):
    assert zf_prime(p) and chain(p) == n, p
print("page chains ok: 11 (1), 1433 (2), 122527 (3), 112222268681 (4)")
# n = 4 minimal: over all multisets of 1..9 with <= 12 digits, the image must be a prime with chain >= 3
qual = []
for size in range(1, 13):
    for ms in itertools.combinations_with_replacement(range(1, 10), size):
        S = sum(POW[d] for d in ms)
        if zf_prime(S) and chain(S) >= 3: qual.append(ms)
assert qual == [(1, 1, 1, 2, 2, 2, 2, 2, 6, 6, 8, 8)], qual
from sympy.utilities.iterables import multiset_permutations
least = min(int("".join(map(str, p))) for p in multiset_permutations(list(qual[0])) if isprime(int("".join(map(str, p)))))
assert least == 112222268681
print("n = 4: the only qualifying multiset with <= 12 digits is 111222226688; its least prime arrangement is 112222268681 (minimal)")

# n = 5 strict
txt = open("writeup/data/p1049_n5_strict.txt").read()
p = int(re.search(r"least p: (\d+)", txt).group(1))
assert len(str(p)) == 321 and '0' not in str(p)
assert pari.isprime(p) == 1
seq = [p]
for _ in range(5): seq.append(f(seq[-1]))
assert all(zf_prime(x) for x in seq[1:]) and not zf_prime(f(seq[-1]))
print("n = 5: p has 321 digits, proven prime (PARI); chain:", " -> ".join(map(str, seq[1:])), "-> ", f(seq[-1]), "(not prime)")
# D(q) for every candidate q by an independent exact method: a DP table T[r] = least number of coins from
# {1, 4, 27, 256, 3125, 46656, 823543} making r, for r < 2 * 8^8; then enumerate c9 and c8 with exact lower
# bounds (every remaining coin is at most 8^8, resp. 7^7), so a state outside the table cannot be optimal.
Ds = {int(m.group(1)): int(m.group(2)) for m in re.finditer(r"q=(\d+) D=(\d+)", txt)}
assert len(Ds) == 5566
R = 2 * POW[8]
T = np.arange(R, dtype=np.int32)              # coins of 1 only
for coin in (4, 27, 256, 3125, 46656, 823543):   # T[r] = min_k T[r - k coin] + k, per residue class by a cumulative minimum
    for s0 in range(coin):
        a = T[s0::coin]; j = np.arange(len(a), dtype=np.int32)
        T[s0::coin] = np.minimum.accumulate(a - j) + j
def D_exact(q):
    best = 10**9
    for c9 in range(q // POW[9], -1, -1):
        r = q - c9 * POW[9]
        if c9 + -(-r // POW[8]) >= best: break          # fewer nines only get worse from here
        for c8 in range(r // POW[8], -1, -1):
            r2 = r - c8 * POW[8]
            if c9 + c8 + -(-r2 // POW[7]) >= best: break
            if r2 < R: best = min(best, c9 + c8 + int(T[r2]))
    return best
for q, D in Ds.items(): assert D_exact(q) == D, (q, D, D_exact(q))
print("D(q) of all 5566 candidates confirmed by the DP-table method; least 321, for q = 112826622281 and 116226228821")
assert f(p) == 112826622281
# the smaller arrangements of the same multiset are composite; the other q's optimal multiset sorts higher
digits = sorted(str(p)); arr = digits[:]; n = 0
while "".join(arr) < str(p):
    n += 1; assert not isprime(int("".join(arr)))
    i = len(arr) - 2
    while arr[i] >= arr[i + 1]: i -= 1
    j = len(arr) - 1
    while arr[j] <= arr[i]: j -= 1
    arr[i], arr[j] = arr[j], arr[i]; arr[i + 1:] = reversed(arr[i + 1:])
assert "".join(arr) == str(p)
other = [m for m in re.finditer(r"q=116226228821 D=321 multisets=\[(.*)\]", txt)][0].group(1)
print(f"{n} smaller arrangements of the multiset are composite; multisets for q = 116226228821: {other}")
# the loose reading (later chain members may contain zeros): data/p1049_loose.txt
lt = open("writeup/data/p1049_loose.txt").read()
assert re.search(r"n=4 loose least p: (\d+)", lt).group(1) == "112222268681"
pl = int(re.findall(r"least p: (\d+)", lt)[-1])
assert len(str(pl)) == 311 and '0' not in str(pl) and pari.isprime(pl) == 1
seq = [pl]
for _ in range(5): seq.append(f(seq[-1]))
assert all(isprime(x) for x in seq[1:]) and not isprime(f(seq[-1]))
zeros = [x for x in seq[1:] if '0' in str(x)]
Dl = {int(m.group(1)): int(m.group(2)) for m in re.finditer(r"q=(\d+) D=(\d+)", lt)}
assert min(Dl.values()) == 311 and D_exact(114693336623) == 311 and f(pl) == 114693336623
for q, D in list(Dl.items())[:300]: assert D_exact(q) == D
print("loose reading: n = 4 least p unchanged; n = 5 least p has 311 digits (proven prime), chain", " -> ".join(map(str, seq[1:])),
      "; members with a zero digit:", zeros, "; D(q) of the first 300 candidates confirmed by the DP method")
print("all checks passed")
