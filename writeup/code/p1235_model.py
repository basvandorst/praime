# Puzzle 1235 (writeup/1235.md): F(n) = 6^n + n^6. How likely is a prime F(n) with n > 1?
# Allowed n: the classes mod 210 not killed by 2, 3, 5, 7: n = 35, 49, 91, 119, 161, 175 (mod 210).
# Local density for a prime q >= 11, conditional on the class c (mod 210):
#   q | F(n) iff 6^n = -n^6 (mod q). Let o = ord_q(6), g = gcd(6, q-1), S = sixth powers (index g), and
#   m = o / gcd(o, (q-1)/g) = the index of H n S in H = <6>. Then 6^s is in -S iff w^s = eps with w = 6^((q-1)/g) of
#   order m and eps = (-1)^((q-1)/g): iff s = 0 (mod m) when eps = 1, iff s = m/2 (mod m) when eps = -1 (m even).
#   For such s there are g values r = n mod q with r^6 = -6^s. The class c fixes n mod d = gcd(o, 210), so among the
#   o/d admissible residues s the condition picks o/lcm(m, d) of them if c is compatible, none otherwise:
#   f_q(c) = g * gcd(m, d) / (q * m) if (eps = 1 and c = 0 mod gcd(m, d)) or (eps = -1 and c = m/2 mod gcd(m, d)),
#   else 0.  Checked by brute force below.
# Model: for n in class c, P(F(n) prime) = (2 * 3/2 * 5/4 * 7/6) * C_c / ln F(n), C_c = prod_{q>=11} (1 - f_q(c))/(1 - 1/q).
# Calibration: survivors of the sieve by 11 <= q <= 1000 among allowed n <= 150000, per class, against the prediction.
# Run: python3 p1235_model.py   (about a minute; needs sympy)
from math import log, gcd, exp
from sympy import primerange, n_order

ALLOWED = [35, 49, 91, 119, 161, 175]
def density(q, c):
    o = n_order(6, q); g = gcd(6, q - 1); m = o // gcd(o, (q - 1) // g); d = gcd(o, 210)
    eps = -1 if ((q - 1) // g) % 2 else 1
    G = gcd(m, d)
    if eps == 1: ok = c % G == 0
    else: ok = m % 2 == 0 and (c - m // 2) % G == 0
    return g * G / (q * m) if ok else 0.0

# brute force check: fraction of n in class c (mod 210) with q | F(n), over a full period lcm(210, q * o)
for q in primerange(11, 50):
    o = n_order(6, q); P = 210 * q * o // gcd(210, q * o)
    for c in ALLOWED:
        cnt = sum(1 for n in range(c, P, 210) if (pow(6, n, q) + pow(n, 6, q)) % q == 0)
        assert abs(cnt / (P // 210) - density(q, c)) < 1e-12, (q, c, cnt / (P // 210), density(q, c))
assert sorted({n % 210 for n in range(210, 420) if all((pow(6, n, q) + pow(n, 6, q)) % q for q in (2, 3, 5, 7))}) == ALLOWED
print('conditional density formula checked by brute force for 11 <= q < 50 and all six classes')

C = {c: 1.0 for c in ALLOWED}; trace = []
for q in primerange(11, 200001):
    for c in ALLOWED: C[c] *= (1 - density(q, c)) / (1 - 1 / q)
    if q in (97, 997, 9973, 99991, 199999): trace.append((q, sum(C.values()) / 6))
print('mean of C_c over the classes, product to Q:', ' '.join('Q=%d: %.4f' % t for t in trace))
print('C_c per class (Q = 2e5):', ' '.join('%d: %.3f' % (c, C[c]) for c in ALLOWED))
small = 2 * 3 / 2 * 5 / 4 * 7 / 6

X = 150000; Q = 1000; qs = list(primerange(11, Q + 1))
tot_s = tot_p = 0
for c in ALLOWED:
    ns = list(range(c, X + 1, 210))
    surv = sum(1 for n in ns if all((pow(6, n, q) + pow(n, 6, q)) % q for q in qs))
    pred = len(ns)
    for q in qs: pred *= 1 - density(q, c)
    tot_s += surv; tot_p += pred
    print('class %3d: %d of %d survive the primes 11..%d, predicted %.1f' % (c, surv, len(ns), Q, pred))
print('calibration: %d survivors, %.1f predicted (ratio %.3f)' % (tot_s, tot_p, tot_s / tot_p))

def expected(X):
    return sum(small * C[c] / (n * log(6) + 6 * log(n)) for c in ALLOWED for n in range(c, X + 1, 210))
for X in (1000, 10000, 150000, 10**6, 10**7):
    print('expected number of primes F(n), 1 < n <= %d: %.3f' % (X, expected(X)))
rate = small * sum(C.values()) / 6 / (log(6) * 35)
print('beyond n = 150000 the expected count grows like %.4f ln(X/150000): one more prime expected near X = 10^%.1f' % (rate, log(150000 * exp(1 / rate), 10)))
print('all checks passed')
