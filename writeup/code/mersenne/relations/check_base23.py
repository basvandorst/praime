"""Independent recheck of the strongest single signal (base 23, longest run of equal
adjacent digits), plain Python, with a third null:
  (c) 'same size' null: each real term a_k (k >= 2) is replaced by a random prime in
      [a_k / sqrt2, a_k * sqrt2] (next prime after a uniform random integer), 5000 samples.
Also the same procedure for p mod 8 == 1 and the binary palindromes."""
import json, math, random
import gmpy2

A = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
     4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049,
     216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
     20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609,
     57885161, 74207281, 77232917, 82589933, 136279841]


def longest_run(x, b):
    d = []
    while x:
        x, r = divmod(x, b)
        d.append(r)
    best = cur = 1
    for i in range(1, len(d)):
        cur = cur + 1 if d[i] == d[i - 1] else 1
        best = max(best, cur)
    return best, len(d)


def run23(seq):
    v = [longest_run(x, 23)[0] for x in seq if x >= 23]
    return sum(v) / len(v), sum(1 for t in v if t >= 2)


def mod8(seq):
    return sum(1 for x in seq if x % 8 == 1)


def binpal(seq):
    return sum(1 for x in seq if x >= 2 and bin(x)[2:] == bin(x)[2:][::-1])


real = (run23(A), run23(A[25:]), mod8(A), binpal(A))
print('real: run23 (mean, #terms with run>=2) =', real[0], ' second half:', real[1],
      ' #p==1 mod 8:', real[2], ' binary palindromes:', real[3])

# same nulls as the main run, recomputed with this code
D = json.load(open('nulls.json'))
for fam in ('a', 'b'):
    vals = [run23(s)[0] for s in D[fam]]
    p = (1 + sum(v >= real[0][0] for v in vals)) / (len(vals) + 1)
    print('null (%s) recomputed: run23 p = %.4f' % (fam, p))

random.seed(2026)
S = 5000
ge = ge2 = ge8 = gepal = 0
cnt2 = []
for _ in range(S):
    seq = [2]
    for x in A[1:]:
        lo, hi = int(x / math.sqrt(2)), int(x * math.sqrt(2))
        while True:
            q = int(gmpy2.next_prime(random.randint(lo, hi)))
            if q <= hi and q > 2:
                break
        seq.append(q)
    r = run23(seq)
    r2 = run23(seq[25:])
    ge += r[0] >= real[0][0]
    ge2 += r2[0] >= real[1][0]
    ge8 += mod8(seq) >= real[2]
    gepal += binpal(seq) >= real[3]
    cnt2.append(r[1])
print('null (c) same-size primes, %d samples:' % S)
print('  run23 full p = %.4f; second half p = %.4f; mean #terms with run>=2 = %.2f (real %d)'
      % ((1 + ge) / (S + 1), (1 + ge2) / (S + 1), sum(cnt2) / S, real[0][1]))
print('  p==1 mod 8 count p = %.4f; binary palindromes p = %.4f' % ((1 + ge8) / (S + 1), (1 + gepal) / (S + 1)))
