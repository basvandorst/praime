"""Puzzle 312, independent check of Q2 with plain integers and sympy.isprime only.
Rubin's rule, reading 'not appeared yet in either sequence': recompute the terms one by one and
show that the sum at the 29th term repeats an earlier sum (463).  Then the same for the reading
'not appeared in the term sequence' (661 at the 49th term)."""
from sympy import isprime
def walk(either):
    P = [3, 5]; sums = []
    while True:
        a = P[-2] + P[-1]
        q = 3
        while q in P or (either and q in sums) or not isprime(q) or not isprime(a + q):
            q += 2
        P.append(q); sums.append(a + q)
        if sums.count(a + q) > 1 or (a + q) in P:
            return P, sums
P, S = walk(True)
assert len(P) == 29 and P[-1] == 157 and S[-1] == 463 and S.index(463) == 24
assert S[24] == P[24] + P[25] + P[26] == 151 + 173 + 139
print("either-sequence reading: 29th term 157 gives the sum 463 = 139 + 167 + 157, already S(25) = 151 + 173 + 139")
P, S = walk(False)
assert len(P) == 49 and P[-1] == 193 and S[-1] == 661 and S.count(661) == 2
print("term-only reading: 49th term 193 gives the sum 661, already used as a sum")
