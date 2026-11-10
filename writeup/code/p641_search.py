# Puzzle 641: chains of semiprimes s -> 2s + 1.  The terms are N_i = 2^i K - 1 with K odd.  K = m * W with
# W = 3*5*...*Y, so that no N_i has a prime factor <= Y.  For each odd m in [M0, M1] the number of semiprimes
# N_1, N_2, ... in a row is counted (PARI bigomega); chains of at least RMIN are printed.
# usage: p641_search.py Y M0 M1 RMIN     e.g.  (ulimit -v 1500000; p641_search.py 13 1 100000000 15)   25 minutes
import sys, cypari2
pari = cypari2.Pari(); pari.allocatemem(200 * 10 ** 6)
Y, M0, M1, RMIN = map(int, sys.argv[1:5])
W = int(pari(f"prod(i=2, primepi({Y}), prime(i))"))
gp = """
chain(W, m0, m1, rmin) = { my(K, i, cnt = vector(40), out = List());
  forstep(m = m0, m1, 2, K = m * W; i = 1;
    while(bigomega(2^i * K - 1) == 2, i++);
    cnt[i]++;
    if(i - 1 >= rmin, listput(out, [m, i - 1])));
  [cnt, Vec(out)] }
"""
pari(gp)
cnt, out = pari(f"chain({W}, {M0 | 1}, {M1}, {RMIN})")
cnt = [int(x) for x in cnt]
tot = sum(cnt); tail = tot
print(f"y = {Y}, W = {W}, odd m in [{M0}, {M1}]: {tot} values")
for r in range(0, 30):
    if tail == 0: break
    print(f"  run >= {r}: {tail}"); tail -= cnt[r]
for m, r in out: print("LONG", int(r), "m", int(m), "K", int(m) * W)
