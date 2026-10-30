"""Puzzle 1238: checks for results/1238.md. Reversed digits come from plain long division here."""
import re
from sympy import n_order, primerange

def period(p):
    """Digits of one period of 1/p by long division."""
    r, out = 1, []
    while True:
        out.append(10 * r // p); r = 10 * r % p
        if r == 1: return out

def reverse_rule_ok(p, N, J, terms):
    """Rotondo's process: start with the first N reversed digits, a(n) = a(n-J) + a(n-N) + carry."""
    rev = period(p)[::-1]
    want = [rev[k % len(rev)] for k in range(terms)]
    d, c = want[:N], 0
    for n in range(N, terms):
        v = d[n - J] + d[n - N] + c
        d.append(v % 10); c = v // 10
    return d == want

# 1. Rotondo's example: 1/7 reversed is 7,5,8,2,4,1 and a(n) = a(n-2) + a(n-3) works from 7,5,8.
assert period(7)[::-1] == [7, 5, 8, 2, 4, 1] and reverse_rule_ok(7, 3, 2, 60)

# 2. Criterion: the rule works iff p | 10^N + 10^J - 1 (all p < 400, N <= 30, over 3 periods + N).
for p in primerange(7, 400):
    T = n_order(10, p)
    for N in range(2, 31):
        for J in range(1, N):
            assert reverse_rule_ok(p, N, J, 3 * T + N) == ((10**N + 10**J - 1) % p == 0), (p, N, J)

# 3. The table: every listed rule works; 'none' exactly when no u, v in <10> mod p have u + v = 1.
for line in open('results/data/1238_reverse_rules.txt'):
    if line.startswith('#'): continue
    t = line.split(); p, T = int(t[0]), int(t[1])
    H = {pow(10, k, p) for k in range(T)}
    solvable = any((1 - u) % p in H for u in H)
    if t[2] == 'none':
        assert not solvable, p
    else:
        N, J = int(t[2]), int(t[3])
        assert solvable and reverse_rule_ok(p, N, J, 3 * T + N) and t[4] == ''.join(map(str, period(p)[::-1][:N]))

# 4. Duality with the forward rules on page 1236: every forward rule with lags {J, N}
#    satisfies 10^-J + 10^-N = 1 (mod p), so lags -J, -N (mod T) give a reverse rule.
gk = re.findall(r"For 1/(\d+), using.*?a\(n\) = a\(n-(\d+)\)\s*\+ a\(n-(\d+)\)",
                open("references/Pu 1236 GK.txt.txt").read(), re.S)
for p, j, n in gk:
    p, j, n = int(p), int(j), int(n)
    assert (pow(10, -j, p) + pow(10, -n, p) - 1) % p == 0
    T = n_order(10, p); J2, N2 = sorted(((-j) % T or T, (-n) % T or T))
    if J2 == N2: N2 += T
    assert (10**N2 + 10**J2 - 1) % p == 0 and reverse_rule_ok(p, N2, J2, 3 * T + N2)
print("all checks pass")
