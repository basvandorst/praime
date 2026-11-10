# Puzzle 1198, targeted family: m = s * q, s built from 2..13 with 31-smooth sigma of each prime power, q a prime with q + 1 31-smooth
# (or q = 1).  For each m the last prime p is solved as in p1198.py.  usage: p1198b.py PARTMAX QMAX G part nparts
import itertools, math, sys, time, cypari2
from sympy import primerange, divisor_sigma
PARTMAX, QMAX, G = int(float(sys.argv[1])), int(float(sys.argv[2])), int(float(sys.argv[3])); part, nparts = int(sys.argv[4]), int(sys.argv[5])
Y = 31
def smooth(n):
    for p in primerange(2, Y + 1):
        while n % p == 0: n //= p
    return n == 1
small_p = [2, 3, 5, 7, 11, 13]; exps = {2: 14, 3: 8, 5: 5, 7: 4, 11: 3, 13: 3}
comp = {p: [e for e in range(0, exps[p] + 1) if e == 0 or smooth(int(divisor_sigma(p ** e)))] for p in small_p}
parts = sorted(x for x in (math.prod(p ** e for p, e in zip(small_p, es)) for es in itertools.product(*[comp[p] for p in small_p])) if x <= PARTMAX)
big = [1] + [q for q in primerange(17, QMAX) if smooth(q + 1)]
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
pari("ssd(n) = { my(f = factor(n)); prod(i = 1, matsize(f)[1], sum(j = 0, f[i,2], sigma(f[i,1]^j))) }")
pari("""scanv(v, G) = { my(out = List(), m, A, B, P, c, num, den, r, p, g, np, lst, nl, x);
  for(t = 1, #v, m = v[t]; A = sigma(m); B = ssd(m); P = factor(A)[,1]; np = #P;
    lst = List([1]);
    for(i = 1, np, nl = #lst; for(j = 1, nl, x = lst[j] * P[i]; while(x <= G, listput(lst, x); x *= P[i])));
    for(k = 1, #lst, g = lst[k]; if(g == 1, next); c = sigma(A * g);
      if(c == B * (g + 1) && isprime(g - 1) && m % (g - 1), listput(out, [m, g - 1]));
      if(c > B && c < B * g, num = c - B; den = B * g - c;
        if(num % den == 0, r = num / den;
          if(r > 1 && gcd(r, A) == 1 && isprime(r), p = g * r - 1; if(isprime(p) && m % p, listput(out, [m, p])))))));
  Vec(out) }""")
print(f"# parts {len(parts)}, primes {len(big)}, G {G}", flush=True)
t0 = time.time(); seen = set()
for i, q in enumerate(big):
    if i % nparts != part: continue
    v = [s * q for s in parts]
    for m, p in pari(f"scanv({v}, {G})"):
        n = int(m) * int(p)
        if n in seen: continue
        seen.add(n); ok = pari(f"ssd({n}) == sigma(sigma({n}))")
        print("SOLUTION" if ok else "false", n, "= m*p", int(m), int(p), pari(f"factor({n})"), flush=True)
    if i % (50 * nparts) == part: print(f"# q index {i} of {len(big)} (q = {q}), {time.time() - t0:.0f} s", flush=True)
print("# done", flush=True)
