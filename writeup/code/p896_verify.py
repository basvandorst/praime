# Puzzle 896 (writeup/896.md): shortest digit string containing every word of a set as a subsequence.
# Checks the strings on the page (set A in 36 digits, primes < 10^4 in 30, primes < 10^5 in 40) and proves with a SAT
# solver that the primes < 100 need 11 digits and the primes < 1000 need 21: no string of length 10, resp. 20, exists.
# Formula: x[i][d] = position i holds digit d (exactly one); for every prefix w_1..w_j of a word and every i,
# m[j][i] <-> m[j][i-1] or (m[j-1][i-1] and x[i][w_j]) ("w_1..w_j is a subsequence of positions 1..i");
# each word asserts m[|w|][L]. Words share prefix variables.
# Run: python3 p896_verify.py   (about 3 minutes; needs sympy, python-sat)
import itertools
from sympy import primerange
from pysat.solvers import Cadical153

def subseq(w, s):
    it = iter(s); return all(c in it for c in w)

def exists(words, L):
    top = [0]
    def new():
        top[0] += 1; return top[0]
    x = [[new() for d in range(10)] for i in range(L)]
    cls = []
    for i in range(L):
        cls.append(x[i][:]); cls += [[-x[i][a], -x[i][b]] for a, b in itertools.combinations(range(10), 2)]
    T = new(); cls.append([T])
    memo = {(): [T] * (L + 1)}
    def m(pre):
        if pre in memo: return memo[pre]
        prev = m(pre[:-1]); d = pre[-1]
        F = new(); cls.append([-F]); cur = [F]
        for i in range(1, L + 1):
            v = new(); a, b, xi = cur[i - 1], prev[i - 1], x[i - 1][d]
            cls.extend([[-a, v], [-b, -xi, v], [-v, a, b], [-v, a, xi]]); cur.append(v)
        memo[pre] = cur; return cur
    for w in words: cls.append([m(tuple(map(int, w)))[L]])
    s = Cadical153(bootstrap_with=cls)
    if not s.solve(): return None
    mod = set(v for v in s.get_model() if v > 0)
    return ''.join(str(next(d for d in range(10) if x[i][d] in mod)) for i in range(L))

A = '4242064951 5519385407 6886477529 3816749731 4700006089 9650560979 4371455527 5419320281 8761823531 3454277693'.split()
for name, words, s in [('set A', A, '435874506198230864527700564095278391'),
                       ('primes < 10^4', [str(p) for p in primerange(2, 10**4)], '425813679450328167906458231793'),
                       ('primes < 10^5', [str(p) for p in primerange(2, 10**5)], '2145879637814520936740281539610872549371')]:
    assert all(subseq(w, s) for w in words)
    print('%-14s the page string of %d digits contains all %d words' % (name, len(s), len(words)))
for name, P, L in [('B: primes < 100', 100, 11), ('C: primes < 1000', 1000, 21)]:
    words = [str(p) for p in primerange(2, P)]
    s = exists(words, L)
    assert s and all(subseq(w, s) for w in words)
    assert exists(words, L - 1) is None
    print('%-17s length %d works (%s), length %d is impossible: %d is the minimum' % (name, L, s, L - 1, L))
print('all checks passed')
