import json
from sympy import primerange
P = list(primerange(2, 200000)); G = [b - a for a, b in zip(P, P[1:])]
from collections import Counter
C = Counter(G[:5000])
json.dump({"gaps200": G[:200], "max_gap_below_2e5": max(G), "counts_first5000": {str(k): v for k, v in sorted(C.items())}}, open("ref40.json", "w"))
print(G[:30], max(G), C.most_common(5))
