# Puzzle 1015 (and the both-orders variant of puzzle 640): sequences of semiprimes n = p*q -> p&q or q&p (concatenation).
# For every semiprime n in [A, B] the longest sequence starting at n is computed by a tree search (PARI; forfactored
# gives the factorisations of the starts).  usage: p1015_search.py A B RMIN     prints the starts with length >= RMIN
# Without arguments: lengths for a few starts of the page's table.
import sys, time, cypari2
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
pari("cat(p, q) = p * 10^#digits(q) + q")
# longest sequence starting at n with factorisation matrix f (n a semiprime)
pari("dep(f) = { my(p = f[1,1], q = if(matsize(f)[1] == 1, p, f[2,1]), g, a = 0, b = 0); g = factor(cat(p, q)); if(bigomega(g) == 2, a = dep(g)); if(p != q, g = factor(cat(q, p)); if(bigomega(g) == 2, b = dep(g))); 1 + max(a, b) }")
pari("scan(a, b, rmin) = { my(cnt = vector(60), out = List(), c, f); forfactored(N = a, b, f = N[2]; if(bigomega(f) == 2, c = dep(f); cnt[c]++; if(c >= rmin, listput(out, [N[1], c])))); [cnt, Vec(out)] }")
if len(sys.argv) == 1:
    print([int(pari(f"dep(factor({n}))")) for n in (6, 4, 15, 10, 178, 35239, 16064229)])
else:
    a, b, rmin = int(float(sys.argv[1])), int(float(sys.argv[2])), int(sys.argv[3])
    t = time.time(); cnt, out = pari(f"scan({a}, {b}, {rmin})"); cnt = [int(x) for x in cnt]
    tot = sum(cnt); print(f"[{a}, {b}]: {tot} semiprimes, {time.time() - t:.1f} s")
    tail = tot
    for r in range(1, 60):
        if tail == 0: break
        if r >= 14 or r % 4 == 1: print(f"  length >= {r}: {tail}")
        tail -= cnt[r - 1]
    for n, c in out: print("LONG", int(c), int(n))
