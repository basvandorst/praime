"""Candidate elliptic curves for the degree-2 test (no elldata on this machine): enumerate Weierstrass equations
[a1,a2,a3,a4,a6] with a1 in {0,1}, a2 in {-1,0,1}, a3 in {0,1}, |a4|, |a6| <= B; keep SEMISTABLE curves (squarefree
conductor, so every bad prime is multiplicative and the bad local factors are exact) with NMIN <= N <= NMAX; one curve per
isogeny class (isogenous curves have the same L-function; key = N and a_p for p < 200); record the global minimal model,
root number, a_p for p <= 13 and the reduction type at 2 and 3.
usage: deg2_curves.py B NMIN NMAX OUT.json"""
import sys, json, time
import cypari2
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
B, NMIN, NMAX, OUT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
PK = [int(p) for p in pari.primes(46)]                       # primes < 200
seen = {}; t0 = time.time(); ntried = 0
for a1 in (0, 1):
    for a2 in (-1, 0, 1):
        for a3 in (0, 1):
            for a4 in range(-B, B + 1):
                for a6 in range(-B, B + 1):
                    ntried += 1
                    try:
                        E = pari.ellinit([a1, a2, a3, a4, a6])
                    except cypari2.PariError:
                        continue
                    if len(E) == 0 or E[11] == 0:
                        continue
                    gr = pari.ellglobalred(E); N = int(gr[0])
                    if N < NMIN or N > NMAX or not bool(pari.issquarefree(N)):
                        continue
                    key = (N,) + tuple(int(pari.ellap(E, p)) for p in PK)
                    if key in seen:
                        continue
                    Em = pari.ellchangecurve(E, gr[1])
                    seen[key] = [int(Em[i]) for i in range(5)]
rows = []
for key, ainv in sorted(seen.items()):
    E = pari.ellinit(ainv); N = key[0]
    rows.append({'ainv': ainv, 'N': N, 'w': int(pari.ellrootno(E)),
                 'ap': {str(p): int(pari.ellap(E, p)) for p in (2, 3, 5, 7, 11, 13)},
                 'bad': [int(p) for p in pari.factor(N)[0]]})
json.dump({'B': B, 'NMIN': NMIN, 'NMAX': NMAX, 'tried': ntried, 'classes': rows}, open(OUT, 'w'))
print('tried %d equations, %d semistable isogeny classes with %d <= N <= %d, %.0fs' % (ntried, len(rows), NMIN, NMAX, time.time() - t0))
