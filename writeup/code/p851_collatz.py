# Puzzle 851 Q2: H(2^n + k) = e(k) + o(k) + 3*floor((n-e)/2) + F_{(n-e)%2}(o + floor((n-e)/2)),
# F_0(A) = H(3^A + 1), F_1(A) = H(2*3^A + 1); H = standard Collatz steps to 1.
import random
def H(x):
    s = 0
    while x != 1:
        x = x // 2 if x % 2 == 0 else 3 * x + 1; s += 1
    return s
def eo(k):
    e = o = 0
    while k != 1:
        if k % 2: k = 3 * k + 1; o += 1
        else: k //= 2; e += 1
    return e, o
n = 1812
# (1) the identity, checked directly on random k
random.seed(1)
for k in [1, 2, 3, 7, 27, 97, 871] + [random.randrange(1, 2**26) for _ in range(40)]:
    e, o = eo(k); t = n - e
    pred = e + o + 3 * (t // 2) + (H(3**(o + t // 2) + 1) if t % 2 == 0 else H(2 * 3**(o + t // 2) + 1))
    assert pred == H(2**n + k), k
print("identity holds for 47 values of k; H(2^1812 + k) =", H(2**n + 1))
# (2) the window: A + F(A) around A = n/2
for par, name in ((0, "3^A+1"), (1, "2*3^A+1")):
    vals = [(A, A + (H(3**A + 1) if par == 0 else H(2 * 3**A + 1))) for A in range(860, 960)]
    const = [A for A, v in vals]
    print(name, "A + F(A) for A = 860..959:", sorted(set(v for A, v in vals)))

# (3) plateaus and the A-window needed by k < 2^m
F0 = {A: H(3**A + 1) for A in range(820, 1000)}
H3 = {A: H(3**A) for A in range(820, 1000)}
plat = []; cur = None
for A in range(820, 1000):
    v = A + F0[A]
    if cur is None or v != cur[2]: cur = [A, A, v]; plat.append(cur)
    else: cur[1] = A
print("plateaus of A + H(3^A+1):", [(a, b, v) for a, b, v in plat])
breaks = [A for A in range(820, 999) if H3[A] != F0[A]]
print("A where H(3^A) != H(3^A + 1) (plateau breaks expected at A+1):", breaks)
import itertools
lo = hi = None
for k in range(1, 2**18):
    e, o = eo(k); A = o + (n - e) // 2
    lo = A if lo is None else min(lo, A); hi = A if hi is None else max(hi, A)
print(f"k < 2^18 reach A in [{lo}, {hi}]")
