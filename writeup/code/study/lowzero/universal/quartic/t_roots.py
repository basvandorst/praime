import cypari2, time
pari = cypari2.Pari(); pari.allocatemem(300_000_000)
f = pari('x^4 - 2*x^2 - 8*x + 1'); P = pari.primes([2, 100000]); print(len(P), pari.nfdisc(f), pari.poldisc(f))
pari('acnt(f, P) = apply(p -> #polrootsmod(f, p) - 1, P)')
pari('acnt2(f, P) = apply(p -> my(g = f*Mod(1,p)); poldegree(gcd(lift(Mod(x, g)^p) - x, g)) - 1, P)')
pari('acnt3(f, P) = apply(p -> my(v=factormod(f,p,1)[,1]); #select(d->d==1, v) - 1, P)')
for name in ('acnt', 'acnt2'):
    t = time.time(); a = pari(name)(f, P); dt = time.time() - t; print(name, '%.2f us per prime' % (1e6*dt/len(P)), list(a[:12]))
