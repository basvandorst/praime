# Puzzle 640 (rule: next term = p&q with p <= q): chain lengths for all semiprimes in [A, B] (PARI).
# usage: p640_stats.py A B RMIN      e.g.  p640_stats.py 8e9 8.004e9 12     without arguments: the page's starts
import sys, time, cypari2
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
pari("nxt(p, q) = p * 10^#digits(q) + q")
pari("len(n) = { my(f, c = 0); while(1, f = factor(n); if(bigomega(f) != 2, return(c)); c++; n = if(matsize(f)[1] == 1, nxt(f[1,1], f[1,1]), nxt(f[1,1], f[2,1]))); }")
pari("scan(a, b, rmin) = { my(cnt = vector(40), out = List(), c, f); forfactored(N = a, b, f = N[2]; if(bigomega(f) == 2, c = 1 + len(if(matsize(f)[1] == 1, nxt(f[1,1], f[1,1]), nxt(f[1,1], f[2,1]))); cnt[c]++; if(c >= rmin, listput(out, [N[1], c])))); [cnt, Vec(out)] }")
if len(sys.argv) == 1:
    print([int(pari(f"len({n})")) for n in (161, 100462, 68279314, 636542506)])
else:
    a, b, rmin = int(float(sys.argv[1])), int(float(sys.argv[2])), int(sys.argv[3])
    t = time.time(); cnt, out = pari(f"scan({a}, {b}, {rmin})"); cnt = [int(x) for x in cnt]
    tot = sum(cnt); print(f"[{a}, {b}]: {tot} semiprimes, {time.time() - t:.1f} s")
    tail = tot
    for r in range(1, 30):
        if tail == 0: break
        print(f"  length >= {r}: {tail}"); tail -= cnt[r - 1] if r - 1 < len(cnt) else 0
    for n, c in out: print("LONG", int(c), int(n))
