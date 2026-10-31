# Puzzle 691: checks the computed claims of results/691.md (about 20 seconds).
from sympy import primerange
P = [int(x) for x in primerange(2, 1_500_000)]
cur = P[:]
first = {}
for k in range(1, 25):
    cur = [b - a for a, b in zip(cur, cur[1:])]
    assert cur[0] % 2 == 1 and all(v % 2 == 0 for v in cur[1:])           # parity (Vantieghem): only the first term is odd
    N = 50000
    assert sum(abs(v) for v in cur[:N]) <= 2**(k - 1) * P[N + k - 1]        # the bound used in the proof
    seen = {}
    for n, v in enumerate(cur, 1):
        if v in seen: first[k] = (seen[v], n, v); break
        seen[v] = n
assert first[9] == (2, 14, 220) and first[12] == (77, 132, 5848) and first[20] == (279, 875, -468924)
assert all(k in first for k in range(1, 25))                                # two equal elements for every k <= 24 already
print({k: first[k] for k in (1, 2, 9, 20, 24)})
print("all claims verified")
