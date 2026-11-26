import sys, time
sys.path.insert(0, sys.argv[1])
from hurwitz import hurwitz12
import cypari2
pari = cypari2.Pari()
for M in (10**5, 10**6, 4*10**6):
    t0 = time.time(); H = hurwitz12(M); t1 = time.time()
    print(M, "%.1fs" % (t1 - t0), flush=True)
    if M == 10**5:
        ref = [int(12 * pari.qfbhclassno(m)) for m in range(1, M + 1)]
        bad = sum(1 for m in range(1, M + 1) if ref[m - 1] != H[m])
        print("mismatches vs PARI for m <= 1e5:", bad)
