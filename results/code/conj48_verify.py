# Conjecture 48 (Sand).  (1) Degiorgi's question: 103 is not |2^n - 3^m|, nor 2^n + 3^m, nor 2^n 3^m +- 1.
# (2) Sand's corrected version (A, B coprime, rad(AB) = product of ALL primes <= Q, nextprime(Q)^2 > P, P = A + B or |A - B|):
#     representations for 397, 419; none for 433 with Q <= 31 and A <= L  (L = 1e13 here; --full uses 1e16, about 20 minutes).
import sys
import numpy as np
from sympy import primerange, factorint
# (1) proof by congruences, step by step
# 2^n + 3^m = 103: 3^m in {1,3,9,27,81} and 103 - 3^m in {102,100,94,76,22}: no power of 2
assert all(((103 - 3**m) & (102 - 3**m)) != 0 for m in range(5))
# 2^n - 3^m = 103: n >= 7; mod 3 -> 2^n = 1 -> n even; mod 4 -> 3^m = 1 -> m even; (2^(n/2) - 3^(m/2))(2^(n/2) + 3^(m/2)) = 103
assert (103 % 3 == 1) and (-103) % 4 == 1 and (103 + 1) // 2 == 52            # 2^(n/2) = (1 + 103)/2 = 52, not a power of 2
# 3^m - 2^n = 103: mod 3 -> 2^n = 2 -> n odd; n >= 3: mod 8 -> 3^m = 7 (impossible, 3^m mod 8 in {1,3}); n = 1: 3^m = 105
assert (-103) % 3 == 2 and 103 % 8 == 7 and {pow(3, m, 8) for m in range(8)} == {1, 3} and 105 != 3**4 and 105 != 3**5
assert not any(abs(2**n - 3**m) == 103 or 2**n + 3**m == 103 or 2**n * 3**m in (102, 104) for n in range(300) for m in range(200))
# (2)
def rep_ok(P, Q, A, B):
    full = set(primerange(2, Q + 1)); fa, fb = set(factorint(A)), set(factorint(B))
    return (A + B == P or abs(A - B) == P) and not (fa & fb) and (fa | fb) == full
assert rep_ok(397, 23, 8008065, 8008462) and rep_ok(419, 29, 12128061, 12128480)
L = 10**16 if "--full" in sys.argv else 10**13
for Q in (19, 23, 29, 31):
    S = list(primerange(2, Q + 1)); sm = [1]
    for p in S:
        new = []
        for x in sm:
            y = x
            while y <= L: new.append(y); y *= p
        sm = new
    arr = np.array(sorted(sm), dtype=np.int64)
    for P in (433,):
        idx = np.searchsorted(arr, arr + P); ok = (idx < len(arr)) & (arr[np.minimum(idx, len(arr) - 1)] == arr + P)
        cands = [(int(a), int(a) + P) for a in arr[ok]] + [(int(a), P - int(a)) for a in arr[arr < P] if (P - int(a)) in set(arr[arr < P].tolist())]
        assert not any(rep_ok(P, Q, A, B) for A, B in cands), (P, Q)
print("433: no representation for Q in {19, 23, 29, 31} with A <=", L)
print("all claims verified")
