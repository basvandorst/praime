from sympy import isprime
def fibmod(m, n):
    a, b = 0, 1 % n
    for bit in bin(m)[2:]:
        c = a * ((2 * b - a) % n) % n; d = (a * a + b * b) % n
        a, b = (d, (c + d) % n) if bit == '1' else (c, d)
    return a
out = []
for k in range(3, 100001):
    if k % 3 == 0 or k % 5 == 0: continue                      # F_k even or divisible by 5
    r5 = k * pow(3, k - 1, 5) % 5                               # F_k mod 5
    eps = 1 if r5 in (1, 4) else -1
    if (fibmod(k, k) - eps) % k == 0:
        kind = "prime" if isprime(k) else ("2p" if k % 2 == 0 and isprime(k // 2) else "other")
        out.append((k, kind))
from collections import Counter
C=Counter(t for _, t in out); print(C)
import json; json.dump({"counts": dict(C), "other": [k for k, t in out if t == "other"]}, open("indices.json", "w"))
oth = [k for k, t in out if t == "other"]
print("other indices k <= 1e5 with k | F_k - eps:", oth)
print("primes k (not 3,5) missing:", [k for k in range(7, 100001) if isprime(k) and k not in dict(out)][:10])
