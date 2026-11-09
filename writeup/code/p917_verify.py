# Puzzle 917: P prime with 2d digits, A and B its two halves (B may start with 0), Q = P + A B, halved if even.
# Checks (1) the examples on the page, (2) every solution of the constructive search in data/p917_solutions.txt
# (all primes proven: sympy.isprime is deterministic below 2^64, and all numbers here are below 10^18),
# (3) every Q with >= 6 preimages found by the exhaustive C sieve over 10-digit primes (data/p917_10digit.txt),
# (4) the exhaustive 8-digit result independently with a numpy sieve (data/p917_8digit.txt).
# Run from the repo root: .venv/bin/python3 writeup/code/p917_verify.py
import re
import numpy as np
from sympy import isprime, divisors, factorint

def Qof(P, d):
    A, B = divmod(P, 10**d); q = P + A * B
    return q // 2 if q % 2 == 0 else q

def check(Q, Ps, d):
    assert len(set(Ps)) == len(Ps) and isprime(Q)
    for P in Ps:
        assert 10**(2*d - 1) <= P < 10**(2*d) and isprime(P) and Qof(P, d) == Q, (Q, P)
    return sum(1 for P in Ps if P % 10**d >= 10**(d - 1))

# (1) the page
check(139801, [102367, 237179, 261071], 3); check(165523, [197677, 233419, 257287, 285161], 3)
check(36513391, [26083999, 37379539, 49134863, 69410521, 70130413], 4)   # Torasso's a) uses 69410521 -> B = 0521
check(39207523, [45217343, 46736779, 47936359, 51655181, 53574637], 4)
check(44017069, [44459803, 46738837, 49137917, 71812259, 74611799], 4)
print("page examples ok (k = 3, 4, 5)")

# (2) constructive solutions
for line in open("writeup/data/p917_solutions.txt"):
    if line.startswith("#"): continue
    v = list(map(int, line.split())); digits, k, ks, Q, Ps = v[0], v[1], v[2], v[3], v[4:]
    assert digits % 2 == 0 and len(Ps) == k and Q < 2**64 and max(Ps) < 2**64
    strict = check(Q, Ps, digits // 2); assert strict == ks
    # identity: Q + 10^d = (A+1)(B+10^d) or 2Q + 10^d = (A+1)(B+10^d)
    d = digits // 2
    for P in Ps:
        A, B = divmod(P, 10**d); assert (A + 1) * (B + 10**d) in (Q + 10**d, 2*Q + 10**d)
    win = lambda N: sum(1 for f in divisors(N) if 10**d <= f < 2 * 10**d and N // f > 10**(d - 1))
    w1, w2 = win(Q + 10**d), win(2 * Q + 10**d)
    print(f"{digits}-digit primes: k = {k} ({ks} without leading zero in B), Q = {Q}: ok; "
          f"Q + 10^d = {factorint(Q + 10**d)} with {w1} divisors in the window, 2Q + 10^d has {w2}")

# (3) the 10-digit sieve output
lines = open("writeup/data/p917_10digit.txt").read().split("\n")
hist = {}
for l in lines:
    m = re.match(r"k=(\d+) Q=(\d+) P:(.*)", l)
    if not m: continue
    k, Q = int(m.group(1)), int(m.group(2)); Ps = [int(x.split("(")[0]) for x in m.group(3).split()]
    assert len(Ps) == k and k >= 6; check(Q, Ps, 5); hist[k] = hist.get(k, 0) + 1
m = re.search(r"primes P with 10 digits: (\d+)", "\n".join(lines)); assert m.group(1) == "404204977"  # pi(1e10) - pi(1e9)
print("10-digit sieve: all listed solutions verified; number of prime Q by k:", dict(sorted(hist.items())))
assert hist == {6: 544, 7: 117, 8: 17, 9: 6, 10: 1}

# (4) 8-digit, independently
N = 10**8; s = np.ones(N, bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.flatnonzero(s[10**7:]).astype(np.int64) + 10**7
A, B = P // 10**4, P % 10**4; q = P + A * B; q = np.where(q % 2 == 0, q // 2, q)
vals, cnt = np.unique(q, return_counts=True)
big = [(int(v), int(c)) for v, c in zip(vals[cnt >= 5], cnt[cnt >= 5]) if isprime(int(v))]
c_lines = [l for l in open("writeup/data/p917_8digit.txt").read().split("\n") if l.startswith("k=5")]
assert sorted(v for v, c in big) == sorted(int(re.search(r"Q=(\d+)", l).group(1)) for l in c_lines) and all(c == 5 for v, c in big)
assert not any(isprime(int(v)) for v in vals[cnt >= 6])   # composite Q may have more
print("8-digit: numpy recount agrees with the sieve: %d prime Q with 5 preimages, none with 6" % len(big))
print("all checks passed")
