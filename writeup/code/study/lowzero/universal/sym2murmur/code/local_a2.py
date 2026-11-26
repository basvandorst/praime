import sys
sys.path.insert(0, sys.argv[1])
from local_a import alpha
for l, emax in ((2, 7), (3, 6), (5, 4), (7, 3)):
    for e in range(1, emax + 1):
        row = []
        for v in range(0, 5):
            a = alpha(l, e, v)
            row.append("%+9.4f" % a.real + ("" if abs(a.imag) < 1e-6 else "(i%+.3f)" % a.imag))
        print("l=%d e=%d:" % (l, e), " ".join(row), flush=True)
