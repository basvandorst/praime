# Joint survival: primes with ord_p(2) | L0 act together on n mod L0 (exact), others independently.
from sympy import primerange, n_order
import sys
L0 = int(sys.argv[1])
obs = {51: (24, 26), 231: (35, 41), 261: (35, 42), 87: (1, 17), 93: (0, 27), 177: (6, 38), 213: (1, 12), 267: (4, 26)}
prs = [(p, n_order(2, p)) for p in primerange(3, 20000)]
for k, (odd, tot) in obs.items():
    alive = [True] * L0; wo = we = 1.0
    for p, o in prs:
        if k % p == 0: continue
        cls = [n for n in range(o) if (k * pow(2, n, p) + 1) % p == 0]
        if not cls: continue
        if L0 % o == 0:
            for r in cls:
                for n in range(r, L0, o): alive[n] = False
        else:
            L = o if o % 2 == 0 else 2 * o
            ko = sum(1 for n in range(L) if n % o in cls and n % 2); ke = sum(1 for n in range(L) if n % o in cls and n % 2 == 0)
            wo *= 1 - ko / (L / 2); we *= 1 - ke / (L / 2)
    so = sum(alive[n] for n in range(1, L0, 2)) / (L0 / 2); se = sum(alive[n] for n in range(0, L0, 2)) / (L0 / 2)
    po, pe = so * wo, se * we
    print(f"k={k:3d}: surviving share of odd n mod {L0}: {so:.4f}, even: {se:.4f};  predicted odd share {po/(po+pe):.3f}   observed (n<=3000) {odd}/{tot} = {odd/tot:.3f}")
