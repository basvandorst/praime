"""All primes of 6 or more rounds up to 2 * 4.2e9, by the lemma instead of a scan (machine 23).
Lemma: each neighbour of a prime of g + 1 rounds has a prime factor of g rounds or more (otherwise that neighbour is
lit completely by round g - 1 and is larger than sqrt N by itself). So a prime N <= 8.4e9 of 6 or more rounds has
N - 1 = s q with q a prime of 5 or more rounds, q <= 4.2e9, and all those q are known (inherit_deep_4e9.txt).
The rounds are worked out from the top down with PARI (inherit_big.py).   usage: inherit_lemma.py [TOP [LIST]]   LIST holds every prime of 5 or more rounds up to TOP / 2 (lines: rounds prime)"""
from inherit_big import rounds, atmost, pari

import sys
TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 8400000000
LIST = sys.argv[2] if len(sys.argv) > 2 else 'inherit_deep_4e9.txt'
deep5 = [int(l.split()[1]) for l in open(LIST) if int(l.split()[1]) <= TOP // 2]
tried = primes = 0
found = {}
for q in deep5:
    s = 2
    while s * q + 1 <= TOP:
        N = s * q + 1
        tried += 1
        if pari.ispseudoprime(N):
            primes += 1
            if not atmost(N, 5):
                found[N] = rounds(N)
        s += 2
print('candidates N = s q + 1 up to %d with q one of the %d primes of 5 or more rounds up to %d: %d, of which %d prime' % (TOP, len(deep5), TOP // 2, tried, primes))
print('primes of 6 or more rounds up to %d: %s' % (TOP, ', '.join('%d (%d rounds)' % (N, g) for N, g in sorted(found.items()))))
print('of these up to 4200000000: %s' % ', '.join(str(N) for N in sorted(found) if N <= 4200000000))
# one level up: a prime of 7 or more rounds up to 2 * TOP is s q + 1 with q one of the primes just found
tried7 = 0
seven = {}
for q in sorted(found):
    s = 2
    while s * q + 1 <= 2 * TOP:
        N = s * q + 1
        tried7 += 1
        if pari.ispseudoprime(N) and not atmost(N, 6):
            seven[N] = rounds(N)
        s += 2
print('candidates N = s q + 1 up to %d with q one of these %d primes: %d; primes of 7 or more rounds up to %d: %s' % (2 * TOP, len(found), tried7, 2 * TOP, ', '.join(str(N) for N in sorted(seven)) or 'none'))

