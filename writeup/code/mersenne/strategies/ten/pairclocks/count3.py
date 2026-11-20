"""Count prime triples (p, p+a, p+b), 7 <= p <= X, by the number of members with an even clock
(H: order of 2 even; L: Fibonacci rank even), reusing the sieve and clocks of count.py.
Usage: python3 count3.py [X]
"""
import json
import sys
import numpy as np
import count as C
from predict3 import PATTERNS

X = C.X


def main():
    res = {p: {"N": 0, **{"%s%d" % (f, j): 0 for f in "HL" for j in range(4)}} for p in PATTERNS}
    for lo in range(0, X + 1, C.SEG):
        hi = min(lo + C.SEG, X + 1)
        s = C.sieve(lo, hi + C.PAD)
        P = (np.nonzero(s)[0] + lo).astype(np.uint64)
        P = P[P >= 7]
        if len(P) == 0:
            continue
        H = C.clock_h(P)
        L = C.clock_l(P)
        Pi = P.astype(np.int64)
        Pm = Pi[Pi < hi]
        for (a, b) in PATTERNS:
            ok = s[Pm + a - lo] & s[Pm + b - lo]
            p1 = Pm[ok]
            idx = [np.searchsorted(Pi, p1 + o) for o in (0, a, b)]
            nh = sum(H[i].astype(int) for i in idx)
            nl = sum(L[i].astype(int) for i in idx)
            r = res[(a, b)]
            r["N"] += len(p1)
            for j in range(4):
                r["H%d" % j] += int((nh == j).sum())
                r["L%d" % j] += int((nl == j).sum())
    json.dump({"X": X, "res": {"%d,%d" % k: v for k, v in res.items()}}, open("count3_%d.json" % X, "w"))
    for k, v in res.items():
        n = v["N"]
        print("(p, p+%d, p+%d), p <= %d: %d triples; " % (k[0], k[1], X, n) +
              ", ".join("%s %.6f" % (q, v[q] / n) for q in v if q != "N"))


if __name__ == "__main__":
    main()
