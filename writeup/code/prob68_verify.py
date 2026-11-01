# Problem 68 (writeup/prob-68.md): Brier numbers k whose minimal covering sets A (for k*2^n + 1) and B (for k*2^n - 1)
# share s primes, with P = product of A u B. Checks the page's records and our examples from scratch:
#   every n mod L has a prime of A dividing k*2^n + 1 and a prime of B dividing k*2^n - 1 (L = lcm of the orders),
#   A and B are minimal (each member is the only one of its set dividing the term for some n),
#   k is odd and larger than every prime used (so each k*2^n +- 1 is a proper multiple: composite),
#   |A n B| and P as claimed.
# Also checks the lemmas of the write-up: a shared prime has even order m and covers n = a (+1 side) and
# n = a + m/2 (-1 side); with 3 shared every other shared prime has m = 2 (mod 4); for m | 720 these are
# 11, 19, 331 and 18837001 only.
# Run: python3 prob68_verify.py   (a few seconds; needs sympy)
from math import lcm, prod
from sympy import isprime, n_order, factorint, cyclotomic_poly

def check(name, k, A, B, s, P):
    A, B = set(A), set(B)
    U = A | B
    assert all(isprime(p) for p in U) and k % 2 == 1 and k > max(U)
    L = lcm(*[n_order(2, p) for p in U])
    privA, privB = set(), set()
    for n in range(L):
        dA = [p for p in A if (k * pow(2, n, p) + 1) % p == 0]
        dB = [p for p in B if (k * pow(2, n, p) - 1) % p == 0]
        assert dA and dB, (name, n)
        if len(dA) == 1: privA.add(dA[0])
        if len(dB) == 1: privB.add(dB[0])
    assert privA == A and privB == B, (name, 'not minimal', A - privA, B - privB)
    assert len(A & B) == s and prod(U) == P, (name, len(A & B), prod(U))
    for p in A & B:                                   # lemma 1: shared primes
        m = n_order(2, p)
        a = [n for n in range(m) if (k * pow(2, n, p) + 1) % p == 0]
        b = [n for n in range(m) if (k * pow(2, n, p) - 1) % p == 0]
        assert m % 2 == 0 and len(a) == len(b) == 1 and (b[0] - a[0]) % m == m // 2
        if p != 3 and 3 in A & B: assert m % 4 == 2   # lemma 2
    print('%-28s s = %d  P = %d (%d digits)  period %d  k = %d' % (name, s, P, len(str(P)), L, k))

# the page (Cohen & Selfridge; Wesolowski 2017 and 2025)
check('page, s = 1 (Cohen-Selfridge)', 47867742232066880047611079,
      [3, 5, 13, 17, 97, 257, 241], [3, 7, 11, 19, 31, 37, 41, 61, 73, 109, 151, 331],
      1, 33241542480794255062005345795)
check('page, s = 2 (Wesolowski)', 1631535784434767401179758185111,
      [3, 5, 17, 19, 97, 109, 241, 257, 433], [3, 7, 11, 13, 19, 31, 37, 41, 61, 73, 151, 331],
      2, 14393587894183912441848314729235)
check('page, s = 3 (Wesolowski)', 1051314309328263323931731864379271,
      [3, 5, 11, 17, 151, 257, 97, 241, 331, 1321, 61681], [3, 7, 11, 13, 19, 31, 37, 41, 61, 73, 109, 331],
      3, 2708540859502146859799487840589422795)
# ours
check('ours, s = 4', 12015917365937901464995049971524255311,
      [3, 7, 11, 13, 19, 37, 61, 73, 109, 151, 331, 1321], [3, 5, 11, 17, 19, 31, 41, 181, 241, 331, 433, 61681],
      4, 8515242199115959558869800655081023415)
check('ours, s = 5', 119644884270899993874221844629387312672591,
      [3, 7, 11, 13, 19, 31, 37, 73, 109, 151, 331, 631, 18837001],
      [3, 5, 11, 17, 19, 41, 61, 97, 181, 241, 257, 331, 18837001],
      5, 71515612535007008102463292424534086157745)
check('ours, s = 6', 98274302749158142265627065837447682246412365955521715066106565731011641111,
      [3, 7, 11, 13, 19, 29, 41, 43, 71, 73, 109, 113, 127, 151, 181, 281, 331, 631, 122921],
      [3, 5, 11, 17, 19, 31, 37, 43, 61, 241, 281, 331, 433, 1429, 3361, 14449, 29191, 106681],
      6, 625440636681326775927408729196212168924063283133237833042367161042005)
check('ours, s = 7', 133719079635465693698861567099435798450160765208807564104962414331041342418162540491,
      [3, 7, 11, 17, 19, 43, 73, 113, 127, 151, 181, 241, 281, 331, 337, 421, 631, 1429, 5419, 61681, 106681, 152041],
      [3, 5, 11, 13, 19, 29, 31, 37, 41, 43, 71, 281, 331, 1321, 5419, 29191, 54001, 7416361],
      7, 138818782164569451915350678846127915041771678931988869163123014322997473719694072155)
check('ours, s = 8', 531778692215922136035345641357241614049447000219099584352521,
      [3, 7, 11, 19, 29, 31, 37, 43, 61, 71, 73, 109, 113, 211, 281, 331, 337, 421, 1429, 5419, 86171],
      [3, 5, 11, 13, 17, 41, 43, 127, 211, 241, 281, 331, 1321, 3361, 5419, 86171],
      8, 279634624029049948951490809829232068694424680986593419861945)

# lemma 3: primes with ord_p(2) = m, m = 2 (mod 4), m | 720
sh = sorted(p for m in range(2, 721, 4) if 720 % m == 0
            for p in factorint(int(cyclotomic_poly(m, 2))) if n_order(2, p) == m)
assert sh == [3, 11, 19, 331, 18837001], sh
print('shareable primes with m = 2 (mod 4), m | 720:', sh)

# the prime list used by prob68_ilp.py is complete for every order m | 720: Phi_m(2) is the product of the listed
# primes of order m (to their full powers) times at most the intrinsic factor, the largest prime of m.
import os, json
here = os.path.dirname(os.path.abspath(__file__))
pr = json.load(open(os.path.join(here, '..', 'data', 'prob68_primes.json')))
for m in range(2, 721):
    if 720 % m: continue
    N = int(cyclotomic_poly(m, 2))
    for p in pr.get(str(m), []):
        assert n_order(2, p) == m and isprime(p) and N % p == 0
        while N % p == 0: N //= p
    assert N == 1 or N == max(factorint(m)), (m, N)
print('prime list complete for all orders m | 720 (%d primes)' % sum(len(pr.get(str(m), [])) for m in range(2, 721) if 720 % m == 0))
big = [p for v in pr.values() for p in v if p >= 2**64]       # sympy's isprime is a proof only below 2^64
try:
    import cypari2
    pari = cypari2.Pari()
    assert all(pari.isprime(p, 2) == 1 for p in big)            # APRCL proof
    print('listed primes above 2^64 proved prime by APRCL:', big)
except ImportError:
    print('cypari2 not installed: primes above 2^64 checked as BPSW probable primes only:', big)
print('all checks passed')
