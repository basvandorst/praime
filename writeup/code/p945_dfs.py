"""Puzzle 945.  For a prime p with digits d_i (d_0 = units), every nonzero digit gives the modulus
d_i 10^i and the remainder r_i = p mod (d_i 10^i) must be prime.  With H_i = the number formed by
the digits above position i and L_i = p mod 10^i:
    r_i = 10^i * ((10 H_i) mod d_i) + L_i.
(10 H_i) mod d_i depends only on H_i mod m(d_i), m = 1,1,3,2,1,3,7,4,9 for d = 1..9; all m divide 252.
Depth-first search from the units digit upward.  State: (i, L_i, M, h) with h = H_i mod M, where M
is the lcm of the moduli still needed; going up, H_i = 10 H_{i+1} + d_{i+1}, so H_{i+1} is known
mod M/gcd(10, M) and the new modulus is lcm(M/gcd(10,M), m(d_{i+1})); the lifts of h are branched
over.  A node can terminate (H_i = 0, i.e. h = 0 mod M) when p = L_{i+1} is prime.
Zero-free digits only.  Enumerates ALL zero-free solutions by length; the tree dies out.
Output: data/p945_solutions.txt (length count smallest largest) and the full list for the longest."""
import sys, time
from math import gcd
from gmpy2 import mpz, is_prime
MOD = {1: 1, 2: 1, 3: 3, 4: 2, 5: 1, 6: 3, 7: 7, 8: 4, 9: 9}
def lcm(a, b): return a * b // gcd(a, b)
sols = {}          # length -> list of primes
maxlen = int(sys.argv[1]) if len(sys.argv) > 1 else 40
nodes = 0
def dfs(i, L, M, h):
    """digits d_0..d_i chosen, L = value of these digits, H_i = h (mod M)."""
    global nodes
    nodes += 1
    # terminate here: H_i = 0 means the top digit is d_i; p = L
    if h % M == 0 and i >= 1 and is_prime(L):
        sols.setdefault(i + 1, []).append(int(L))
    if i + 1 >= maxlen: return
    # next digit d at position i+1: H_i = 10 H_{i+1} + d  ->  10 H_{i+1} = h - d (mod M)
    g = gcd(10, M); M1 = M // g
    for d in range(1, 10):
        if (h - d) % g: continue
        inv = pow(10 // g, -1, M1) if M1 > 1 else 0
        h1 = ((h - d) // g * inv) % M1 if M1 > 1 else 0      # H_{i+1} mod M1
        m = MOD[d]; M2 = lcm(M1, m)
        pw = mpz(10) ** (i + 1)
        for lift in range(h1, M2, M1):                       # lifts of H_{i+1} to Z/M2
            r = pw * ((10 * lift) % d) + L                   # remainder at position i+1
            if not is_prime(r): continue
            dfs(i + 1, L + d * pw, M2, lift)
t0 = time.time()
for d0 in (3, 7, 9):                                          # units digit; its remainder is (10 H) mod d0
    m = MOD[d0]
    for h in range(m):
        if is_prime((10 * h) % d0):
            dfs(0, mpz(d0), m, h)
print(f"nodes {nodes}, {time.time()-t0:.0f}s")
with open("writeup/data/p945_solutions.txt", "w") as f:
    f.write("# puzzle 945: zero-free primes whose remainders modulo every place value d*10^i are prime\n# length count smallest largest\n")
    for n in sorted(sols):
        s = sorted(set(sols[n])); f.write(f"{n} {len(s)} {s[0]} {s[-1]}\n")
        print(n, len(s), s[0], s[-1])
    top = max(sols); f.write("# all solutions of the two largest lengths\n")
    for n in sorted(sols)[-2:]:
        for p in sorted(set(sols[n])): f.write(f"{n} {p}\n")
