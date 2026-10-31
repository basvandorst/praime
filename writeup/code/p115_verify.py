"""Puzzle 115: primes whose square uses at most two different digits.
 (1) exact count for squares below 10^14 (roots < 10^7): all non-trivial (A not divisible by 10) and prime roots;
 (2) the digit-DP model (p115_model.py) for the same lengths, then the expected number above Mishima's bound 10^48;
 (3) the digit pairs possible for p^2 (p prime > 11) at lengths 15..40, from N mod 83160 alone.
Run from writeup/code; about 15 minutes."""
import math, itertools, numpy as np
from sympy import isprime
import p115_model as m
known = [2, 3, 5, 7, 11, 109, 173, 81619]
found = [A for A in range(2, 10**7) if len(set(str(A * A))) <= 2 and isprime(A)]
assert found == known, found
nontriv = {}
for A in range(1, 10**7):
    if A % 10 and len(set(str(A * A))) <= 2: nontriv[len(str(A * A))] = nontriv.get(len(str(A * A)), 0) + 1
u = np.arange(m.M); keep = (u % 10) != 0
rho_all = m.rho.copy(); m.rho = np.bincount((u[keep] * u[keep]) % m.M, minlength=m.M).astype(float)
model_nt = sum(m.expected(k, False) for k in range(5, 15)); m.rho = rho_all
model_p = sum(m.expected(k, True) for k in range(5, 15))
print(f"lengths 5..14: non-trivial squares exact {sum(nontriv.get(k, 0) for k in range(5, 15))}, model {model_nt:.2f}; "
      f"prime squares exact 3 (11881, 29929, 6661661161), model {model_p:.2f}")
E = {k: m.expected(k, True) for k in range(15, 61)}
r = (E[59] + E[60]) / (E[57] + E[58])
print(f"expected prime squares of 15..48 digits (searched, none): {sum(E[k] for k in range(15, 49)):.4f}")
print(f"expected above 10^48: {sum(E[k] for k in range(49, 61)) + (E[59] + E[60]) * r / (1 - r):.2e}")
pairs = set()
for k in range(15, 41):
    for a, b in itertools.combinations(range(10), 2):
        tail = m.dp_counts((a, b), k - 1)
        if any(np.roll(tail, lead * m.pow10[k - 1] % m.M)[m.Rmask].sum() > 0 for lead in (a, b) if lead): pairs.add(f"{a}{b}")
assert sorted(pairs) == ['01', '12', '14', '16', '18', '29', '49', '89'], pairs
print("possible digit pairs for p^2 at 15..40 digits:", " ".join(sorted(pairs)))
