import sys
import numpy as np

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000


def ghaffor_union_size(n):
    top = 6 * n * n + 10 * n + 1
    hit = np.zeros(top + 1, dtype=bool)
    for k in range(1, n + 1):
        q1, q2 = 6 * k - 1, 6 * k + 1
        for first, q in ((5 * k - 3, q1), (7 * k - 3, q1), (5 * k - 1, q2), (7 * k - 1, q2)):
            hit[first::q] = True
    hit[0] = False
    return int(hit.sum())


def ghaffor_count(n):
    return 6 * n * n + 10 * n + 4 - ghaffor_union_size(n)


def twin_counts_upto(limit):
    s = np.ones(limit + 3, dtype=bool)
    s[:2] = False
    for i in range(2, int((limit + 2) ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    twin_lesser = s[:-2] & s[2:]
    return np.cumsum(twin_lesser)


print("Ghaffor's worked examples: |U_1|, |U_2|, |U_3| =", [ghaffor_union_size(n) for n in (1, 2, 3)], "(paper: 10, 29, 63)")
print("G_1, G_2, G_3 =", [ghaffor_count(n) for n in (1, 2, 3)], "(paper: 10, 19, 25)")

limit = 36 * NMAX * NMAX + 60 * NMAX + 21
cum = twin_counts_upto(limit)
bad = []
for n in range(1, NMAX + 1):
    N = 36 * n * n + 60 * n + 21
    if ghaffor_count(n) != int(cum[N]):
        bad.append(n)
print(f"G_n equals the twin count up to 36n^2+60n+21 for all n <= {NMAX}: {not bad}" + (f"  failures {bad[:5]}" if bad else ""))

ok = True
for n in range(1, 61):
    top = 6 * n * n + 10 * n + 1
    hit = np.zeros(top + 1, dtype=bool)
    for k in range(1, n + 1):
        q1, q2 = 6 * k - 1, 6 * k + 1
        for first, q in ((5 * k - 3, q1), (7 * k - 3, q1), (5 * k - 1, q2), (7 * k - 1, q2)):
            hit[first::q] = True
    for v in range(1, top + 1):
        m = v + 2
        a, b = 6 * m - 1, 6 * m + 1
        composite = any(a % d == 0 for d in range(2, int(a**0.5) + 1)) or any(b % d == 0 for d in range(2, int(b**0.5) + 1))
        if hit[v] != composite:
            ok = False
print(f"U_n (index v) = {{v : 6(v+2)-1 or 6(v+2)+1 composite}} for all n <= 60: {ok}")

forms = set()
for k in range(1, 40):
    for b in range(0, 40):
        forms.add((5 * k - 3 + b * (6 * k - 1)) + 2)
        forms.add((7 * k - 3 + b * (6 * k - 1)) + 2)
        forms.add((5 * k - 1 + b * (6 * k + 1)) + 2)
        forms.add((7 * k - 1 + b * (6 * k + 1)) + 2)
perry = set()
for x in range(1, 40):
    for y in range(1, 41):
        perry.update({6 * x * y + x + y, 6 * x * y - x - y, 6 * x * y + x - y, 6 * x * y - x + y})
M = 6 * 39 * 39 - 39 - 40
print(f"Ghaffor's terms (as m = v + 2) and Perry's 6xy ± x ± y agree below {M}: {set(f for f in forms if f < M) == set(p for p in perry if p < M)}")
