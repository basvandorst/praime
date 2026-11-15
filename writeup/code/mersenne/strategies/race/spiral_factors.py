"""Local factors (1 - w(p)/p)/(1 - 1/p) for p <= 13 of selected rays, by brute-force root counting,
and the discriminant of 2f. Also: parity rule for odd P, the s = -1 factorisation."""
import math
rays = [(38, 38, -2), (34, 10, 6), (34, 24, 6), (18, 18, -2), (16, 2, -2), (10, 10, 0), (4, 2, -2), (6, 6, 0), (12, 12, 0), (8, 4, 6)]
def f(P, j, s, m): return P*m*(m-1)//2 + j*m + 1 + s
for P, j, s in rays:
    D = (2*j - P)**2 - 8*P*(1+s)
    fac = []
    prod = 1.0
    for p in (2, 3, 5, 7, 11, 13):
        per = 4*p if p == 2 else p
        w = sum(1 for m in range(per) if f(P, j, s, m) % p == 0) * p / per
        x = (1 - w/p)/(1 - 1/p)
        prod *= x
        fac.append(f"{p}:{x:.3f}")
    print(f"P {P:2d} j {j:2d} s {s:+d}  disc(2f) = {D:6d}  local factors {' '.join(fac)}  product(p<=13) = {prod:.3f}")
# odd P: exactly half the values are even over every period of 4
ok = all(sum(f(P, j, s, m) % 2 for m in range(4)) == 2 for P in range(3, 41, 2) for j in range(1, P+1) for s in range(-6, 7))
print("odd P: f(m) odd for exactly 2 of every 4 consecutive m, all rays:", ok)
# s = -1: f(m) = m * (P(m-1)/2 + j)  -> 2f = m (P(m-1) + 2j)
ok = all(2*f(P, j, -1, m) == m*(P*(m-1) + 2*j) for P in range(3, 41) for j in range(1, P+1) for m in range(1, 200))
print("s = -1: 2 f(m) = m (P(m-1) + 2j) for all P, j, m < 200:", ok)
