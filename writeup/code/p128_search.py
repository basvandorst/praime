# Puzzle 128, the two open parts.  Q[n] = p_1^2 + ... + p_n^2 (int64, exact below 2^62).
# (a) Is Q[n] (start 2) or Q[n] - 4 (start 3) a square?  Exact test for all p_n <= PA, plus the heuristic
#     expectation sum 1/(2 sqrt(N)) over the tested sums and its tail beyond a bound B: sqrt(3)/sqrt(B ln B).
# (c) Two adjacent windows of consecutive primes with equal sums of squares: Q[j] - Q[i] = Q[l] - Q[j], i.e.
#     2 Q[j] - Q[i] is again a prefix sum.  All i < j <= J are tested by binary search in Q; the expectation is
#     sum over the pairs of 1/p_l^2 (the prefix sums near Q[l] are p_l^2 apart).
#     The same test on first powers (J1 primes) must find solutions: a check of the method.
# Usage: p128_search.py [J]     default J = 100000 (110 s, 90 MB); run under ulimit -v 1500000
import sys, math
import numpy as np
J = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
def primes_upto(n):
    sv = np.ones(n + 1, dtype=bool); sv[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if sv[i]: sv[i * i::i] = False
    return np.flatnonzero(sv)

# ---- (a)
PA = 5_000_000
P = primes_upto(PA); Q = np.cumsum(P ** 2)
assert int(P[-1]) ** 2 * len(P) < 2 ** 63                             # no overflow in the prefix sums
for name, N in (("p_1 = 2", Q[1:]), ("p_1 = 3", Q[2:] - 4)):          # at least two primes
    r = np.rint(np.sqrt(N.astype(np.float64))).astype(np.int64)
    sq = np.zeros(len(N), dtype=bool)
    for d in (-1, 0, 1): sq |= (r + d) * (r + d) == N
    e = float((0.5 / np.sqrt(N.astype(np.float64))).sum())
    print(f"(a) {name}: {len(N)} sums with last prime <= {P[-1]}: {int(sq.sum())} squares; expected {e:.3f}, of which {float((0.5 / np.sqrt(N[:10].astype(np.float64))).sum()):.3f} from the first 10 sums")
B = 1693930336951                                                     # Resta's bound on the last prime (page)
print(f"(a) heuristic tail beyond last prime {B}: {math.sqrt(3) / math.sqrt(B * math.log(B)):.1e} per start")
print(f"(a) heuristic tail beyond last prime {P[-1]}: {math.sqrt(3) / math.sqrt(int(P[-1]) * math.log(int(P[-1]))):.1e} per start")

# ---- (c)
def adjacent(k, J):
    n_est = int(2.6 * J * (math.log(J) + math.log(math.log(J)))) + 1000   # second window ends below 2 p_J
    P = primes_upto(n_est)
    assert P[J] * 2 < P[-1]
    Q = np.concatenate(([0], np.cumsum(P ** k)))
    assert Q[-1] < 2 ** 62
    sols = []; exp = 0.0; half = 0.0
    for j in range(1, J + 1):                                         # windows P[i:j] and P[j:l]
        T = 2 * Q[j] - Q[:j]                                          # the prefix sum Q[l] that is needed
        assert T[0] <= Q[-1]
        pos = np.searchsorted(Q, T)
        for i in np.flatnonzero(Q[pos] == T):
            l = int(pos[i]); sols.append((int(P[i]), int(P[j - 1]), j - int(i), int(P[j]), int(P[l - 1]), l - j))
        exp += float((1.0 / P[np.minimum(pos, len(P) - 1)].astype(np.float64) ** k).sum())
        if j == J // 2: half = exp
    return sols, exp, half, int(P[J - 1])
sols, exp, half, last = adjacent(1, 3000)
print(f"(c) first powers, first window ending at p_j <= {last}: {len(sols)} solutions, expected {exp:.1f}; e.g. {sols[0]}")
for a, b, La, c, d, Lb in sols[:3]:
    assert La > 0 and Lb > 0
sols, exp, half, last = adjacent(2, J)
print(f"(c) squares, first window ending at p_j <= {last} (j <= {J}): {len(sols)} solutions, expected {exp:.3f}")
for s in sols: print("   ", s)
# sum over j of j/p_j^2 behaves like the sum of 1/(j ln^2 j), whose tail is 1/ln J: fit the constant on (J/2, J]
c = (exp - half) / (1 / math.log(J // 2) - 1 / math.log(J))
print(f"(c) expected with j in ({J // 2}, {J}]: {exp - half:.4f} = {c:.3f} (1/ln {J // 2} - 1/ln {J}); tail beyond J: {c:.3f}/ln J = {c / math.log(J):.3f}")
