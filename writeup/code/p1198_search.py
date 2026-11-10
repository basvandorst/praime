# Puzzle 1198: ssd(n) = sigma(sigma(n)), ssd(n) = sum of sigma(d) over d | n.  Search for n = m * p, p a prime not dividing m.
# sigma(n) = A (p + 1), ssd(n) = B (p + 2) with A = sigma(m), B = ssd(m).  Write p + 1 = g h, g built from the primes of A,
# gcd(h, A) = 1.  Then sigma(A g) sigma(h) = B (g h + 1).  h = 1 is tested; h = r prime gives r = (c - B)/(B g - c), c = sigma(A g).
# usage: p1198_search.py M0 M1 G    every m in [M0, M1], every g <= G
import sys, time, cypari2
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
pari("ssd(n) = { my(f = factor(n)); prod(i = 1, matsize(f)[1], sum(j = 0, f[i,2], sigma(f[i,1]^j))) }")
pari("""scan(m0, m1, G) = { my(out = List(), A, B, P, c, num, den, r, p, gs, g, np, lst, nl, x);
  for(m = m0, m1, A = sigma(m); B = ssd(m); P = factor(A)[,1]; np = #P;
    lst = List([1]); 
    for(i = 1, np, nl = #lst; for(j = 1, nl, x = lst[j] * P[i]; while(x <= G, listput(lst, x); x *= P[i])));
    for(k = 1, #lst, g = lst[k]; if(g == 1, next); c = sigma(A * g);
      if(c == B * (g + 1) && isprime(g - 1) && m % (g - 1), listput(out, [m, g - 1]));
      if(c > B && c < B * g, num = c - B; den = B * g - c;
        if(num % den == 0, r = num / den;
          if(r > 1 && gcd(r, A) == 1 && isprime(r), p = g * r - 1; if(isprime(p) && m % p, listput(out, [m, p])))))));
  Vec(out) }""")
m0, m1, G = int(float(sys.argv[1])), int(float(sys.argv[2])), int(float(sys.argv[3]))
STEP = 250000
for a in range(m0, m1 + 1, STEP):
    b = min(a + STEP - 1, m1)
    t = time.time(); out = pari(f"scan({a}, {b}, {G})")
    for m, p in out:
        n = int(m) * int(p); ok = pari(f"ssd({n}) == sigma(sigma({n}))")
        print("SOLUTION" if ok else "false", n, "= m*p", int(m), int(p), pari(f"factor({n})"), flush=True)
    print(f"# m in [{a}, {b}], g <= {G}: {time.time() - t:.1f} s", flush=True)
