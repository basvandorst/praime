# Puzzle 845: runs of consecutive integers with an even number of prime factors (Omega, with multiplicity).
# Checks data/p845_runs.txt (output of p845_runs.c): for every reported first occurrence of a run of exactly L,
# recompute Omega of the L members and of the two neighbours with sympy.factorint; compare with the OEIS terms
# A275508(1..30) and the record list A066793 where they are below the search bound.
# Run from the repo root: .venv/bin/python3 writeup/code/p845_verify.py
import re
from sympy import factorint
def even(n): return sum(factorint(n).values()) % 2 == 0
txt = open("writeup/data/p845_runs.txt").read()
N = int(re.search(r"N = (\d+)", txt).group(1))
first = {int(L): int(s) for L, s in re.findall(r"exactly (\d+): (\d+)", txt)}
for L, s in sorted(first.items()):
    assert all(even(s + i) for i in range(L)), (L, s)
    assert (s == 1 or not even(s - 1)) and not even(s + L), (L, s)
print(f"{len(first)} first occurrences below {N} verified by factorisation (exact run lengths)")
A275508 = [1,9,14,33,54,201,140,2170,213,4529,10403,1934,35811,162144,38405,414225,200938,389409,1792209,5606457,
           18493931,30947532,61190563,8405437,203631499,577699919,158280474,68780189,4082994208,3944563444]
assert all(first[L] == A275508[L - 1] for L in range(1, 31)), "mismatch with A275508"
print("A275508(1..30) reproduced")
A066793 = [1,9,14,33,54,140,213,1934,35811,38405,200938,389409,1792209,5606457,8405437,68780189,880346227,
           85910903650,87635944146,284340107954,2195312273836,2864568747364,3155591650029]
recs = [(int(s), int(L)) for s, L in re.findall(r"\((\d+), (\d+)\)", txt)]
assert [s for s, L in recs] == [a for a in A066793 if a < N], (recs, A066793)
print("record runs below the bound agree with A066793:", recs)
new = {L: s for L, s in first.items() if L > 30}
print("new first occurrences (L > 30):", new)
print("all checks passed")
