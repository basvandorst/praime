# Puzzle 945: zero-free primes p such that p mod (d 10^i) is prime for every digit d (at position i) of p.
# (1) Direct definition check of every solution listed in data/p945_solutions.txt (the smallest and largest
#     of each length, and all solutions of the two largest lengths); primality with sympy (proof below 2^64,
#     the 20-23 digit numbers are checked with PARI isprime as well).
# (2) Brute force over all zero-free primes below 10^8 with the direct definition; the counts per length
#     must equal those of the DFS (p945_dfs.py).  (3) The page's values.
# Run from the repo root: .venv/bin/python3 writeup/code/p945_verify.py
from sympy import isprime, primerange
import cypari2
pari = cypari2.Pari()

def ok(p):
    s = str(p)
    if '0' in s or not isprime(p): return False
    n = len(s)
    for i, ch in enumerate(reversed(s)):          # i = position (units = 0), ch = digit
        d = int(ch); r = p % (d * 10**i)
        if not isprime(r): return False
    return True

rows = {}; full = {}
for line in open("writeup/data/p945_solutions.txt"):
    if line.startswith("#"): continue
    v = line.split()
    if len(v) == 4: rows[int(v[0])] = (int(v[1]), int(v[2]), int(v[3]))
    else: full.setdefault(int(v[0]), []).append(int(v[1]))
for n, (cnt, lo, hi) in rows.items():
    assert len(str(lo)) == len(str(hi)) == n and ok(lo) and ok(hi), n
for n, ps in full.items():
    assert len(ps) == rows[n][0] and ps[0] == rows[n][1] and ps[-1] == rows[n][2]
    for p in ps: assert ok(p) and pari.isprime(p) == 1
print("direct check passed for the extremes of every length and all solutions of lengths", sorted(full))

# page values: Izquierdo's list (smallest with 2..11 digits), Resta's smallest 2..16, Cleary 726862627433, Formusa 49963422857
resta = [17, 113, 1223, 12113, 121283, 1214657, 12336353, 121847393, 1215614657, 12184279817, 121681658963,
         1215613379633, 12751425844733, 121947642682373, 1219213679413217]
for n, p in zip(range(2, 17), resta): assert rows[n][1] == p, (n, p)
assert ok(726862627433) and ok(49963422857) and 49963422857 <= rows[11][2] and 726862627433 <= rows[12][2]
print("page values ok: Resta's smallest of each length 2..16 equal ours; Cleary's and Formusa's primes are solutions")

# brute force below 1e8
counts = {}
for p in primerange(2, 10**8):
    if ok(p): counts[len(str(p))] = counts.get(len(str(p)), 0) + 1
print("brute force counts below 1e8:", counts)
assert all(counts.get(n, 0) == rows[n][0] for n in range(2, 9)) and counts.get(1, 0) == 0
print("total zero-free solutions:", sum(c for c, _, _ in rows.values()), "longest:", max(rows), "largest:", rows[max(rows)][2])
print("all checks passed")
