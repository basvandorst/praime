# Puzzle 582 (follow-up question), five prime factors.  With signed c_j, prod (c_j + k i) has imaginary part k * E,
# E = e4 - k^2 e2 + k^4 (e_j = elementary symmetric functions of c_1..c_5); m = a^2 + k^2 needs E = +-1.
# E is linear in each c_j: the four smallest primes are enumerated (all signs), the largest is solved.
# Complete for five-prime products m <= M.   usage: p582_five.py k M      e.g.  p582_five.py 7 1e25
import sys, math, numpy as np, gmpy2
k = int(sys.argv[1]); M = int(float(sys.argv[2])); K2 = k * k
pmin = next(c * c + K2 for c in range(1, 10 ** 5) if gmpy2.is_prime(c * c + K2))
cmax = math.isqrt(math.isqrt(M // pmin ** 3)) + 1
good = np.array([c for c in range(1, cmax + 1) if gmpy2.is_prime(c * c + K2)], dtype=np.int64)
gp = good * good + K2
best = None; ntup = 0
for i1, c1 in enumerate(good.tolist()):
    p1 = c1 * c1 + K2
    if p1 ** 5 > M: break
    for i2 in range(i1, len(good)):
        p2 = int(gp[i2])
        if p1 * p2 ** 4 > M: break
        for i3 in range(i2, len(good)):
            p3 = int(gp[i3])
            if p1 * p2 * p3 ** 3 > M: break
            lim = math.isqrt(M // (p1 * p2 * p3))
            hi = int(np.searchsorted(gp, lim, side="right"))
            if hi <= i3: continue
            c4 = good[i3:hi]; ntup += 8 * len(c4)
            for s2 in (1, -1):
                for s3 in (1, -1):
                    a, b, c = c1, s2 * int(good[i2]), s3 * int(good[i3])
                    f1, f2, f3 = a + b + c, a * b + a * c + b * c, a * b * c
                    for s4 in (1, -1):
                        d = s4 * c4
                        e1 = f1 + d; e2 = f2 + d * f1; e3 = f3 + d * f2; e4 = d * f3
                        A = e4 - K2 * e2 + K2 * K2; B = e3 - K2 * e1
                        for t in (1, -1):
                            nz = B != 0
                            ok = nz & ((t - A) % np.where(nz, B, 1) == 0)
                            for j in np.flatnonzero(ok):
                                x = (t - int(A[j])) // int(B[j])
                                if abs(x) >= int(c4[j]) and gmpy2.is_prime(x * x + K2):
                                    m = p1 * p2 * p3 * (int(c4[j]) ** 2 + K2) * (x * x + K2)
                                    if m <= M and (best is None or m < best[0]): best = (m, (a, b, c, int(d[j]), x))
print("k", k, "M", f"{M:.1e}", "good c up to", cmax, len(good), "tuples", ntup, "best", best, flush=True)
if best:
    m = best[0]; r = math.isqrt(m - K2); assert r * r + K2 == m; print("  a =", r, "primes", [c * c + K2 for c in best[1]])
