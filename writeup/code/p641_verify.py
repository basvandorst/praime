# Puzzle 641: chains of semiprimes s -> 2s + 1.  Independent check with sympy.factorint.
# For each start: every term is a semiprime, (s(1) - 1)/2 is not a semiprime (or not an integer), 2 s(last) + 1 is not.
# Also the page's chains (14, 15, 16 and 17 terms).
# Run from the repo root: .venv/bin/python3 writeup/code/p641_verify.py [start length ...]
import sys
from sympy import factorint
def omega(n): return sum(factorint(n).values())
def chain(s):
    out = []
    while omega(s) == 2: out.append(s); s = 2 * s + 1
    return out
tests = [("page", 2681669, 14), ("page", 121360769, 15), ("page", 969870257, 16), ("page", 9932943134, 17)]
args = [int(x) for x in sys.argv[1:]]
tests += [("ours", args[i], args[i + 1]) for i in range(0, len(args), 2)] or [("ours", 1415485281209, 18)]
for who, s, L in tests:
    c = chain(s)
    before = (s - 1) % 2 == 0 and omega((s - 1) // 2) == 2
    assert len(c) == L and not before, (who, s, len(c), before)
    print(f"{who}: {L} semiprimes from {s}; (s - 1)/2 is not a semiprime; 2*{c[-1]} + 1 has {omega(2 * c[-1] + 1)} prime factors")
    if who == "ours":
        for x in c:
            f = factorint(x); print("   ", x, "=", " * ".join(str(p) for p, e in sorted(f.items()) for _ in range(e)))
print("all checks passed")
