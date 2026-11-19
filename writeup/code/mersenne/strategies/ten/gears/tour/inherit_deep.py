"""Looks for primes that need more rounds than any prime below 4.2e9 (machine 23).
Each neighbour of a prime of g + 1 rounds has a prime factor of g rounds or more. So:
stage 1: q = s r + 1 and s r - 1 for the primes r of 5 or more rounds below 4.2e9 (inherit_deep_4e9.txt) and even s <= S1;
         keep the q that need 6 rounds or more;
stage 2: N = s q + 1 and s q - 1 for these q and even s <= S2; report the N that need 7 rounds or more.
The rounds are worked out from the top down with PARI (inherit_big.py). This finds examples, not the smallest one.
usage: inherit_deep.py S1 S2 [QMAX]      with QMAX: stage 1 takes s up to min(S1, QMAX / r, r) and writes inherit_deep_six_wide.txt
       inherit_deep.py stage2 BOUND [FILE] every candidate s q +- 1 below BOUND for the q of FILE (default inherit_deep_six.txt)"""
import sys
from inherit_big import rounds, atmost, pari, FAC, MEMO

if sys.argv[1] == 'stage2':
    # every candidate s q +- 1 below BOUND, q from inherit_deep_six.txt
    BOUND = int(sys.argv[2])
    six = [int(x) for x in open(sys.argv[3] if len(sys.argv) > 3 else 'inherit_deep_six.txt').read().split()]
    tried = primes = 0
    found = []
    for q in six:
        for s in range(2, BOUND // q + 2, 2):
            for N in (s * q + 1, s * q - 1):
                if N >= BOUND:
                    continue
                tried += 1
                if pari.ispseudoprime(N):
                    primes += 1
                    if not atmost(N, 6):
                        found.append((N, rounds(N), s, q))
        if len(FAC) > 200000:
            FAC.clear()
            MEMO.clear()
    found.sort()
    print('stage 2 below %d: %d primes q of 6 rounds or more, %d candidates s q +- 1, %d primes, %d need 7 rounds or more' % (BOUND, len(six), tried, primes, len(found)))
    for N, g, s, q in found[:12]:
        print('  %d needs %d rounds (%d * %d %s 1)' % (N, g, s, q, '+' if N == s * q + 1 else '-'))
    print('done')
    sys.exit()
S1, S2 = int(sys.argv[1]), int(sys.argv[2])
QMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 0
deep5 = [int(l.split()[1]) for l in open('inherit_deep_4e9.txt')]
print('primes of 5 or more rounds below 4.2e9: %d' % len(deep5))
sys.stdout.flush()
six = set()
tried = primes = 0
for i, r in enumerate(deep5):
    for s in range(2, (min(S1, QMAX // r + 1, r) if QMAX else S1) + 1, 2):
        for q in (s * r + 1, s * r - 1):
            tried += 1
            if pari.ispseudoprime(q):
                primes += 1
                if not atmost(q, 5):
                    six.add(q)
                    g = rounds(q)
                    if g >= 7:
                        print('stage 1 already found %d rounds: %d = %d * %d %s 1' % (g, q, s, r, '+' if q == s * r + 1 else '-'))
                        sys.stdout.flush()
    if len(FAC) > 200000:
        FAC.clear()
        MEMO.clear()
six = sorted(six)
print('stage 1: %d candidates s r +- 1 with s <= %d, %d primes, %d of them need 6 rounds or more; smallest %s' % (tried, S1, primes, len(six), six[:6]))
open('inherit_deep_six%s.txt' % ('_wide' if QMAX else ''), 'w').write('\n'.join(str(q) for q in six) + '\n')
sys.stdout.flush()
tried = primes = 0
found = []
for q in six:
    for s in range(2, S2 + 1, 2):
        for N in (s * q + 1, s * q - 1):
            tried += 1
            if pari.ispseudoprime(N):
                primes += 1
                if not atmost(N, 6):
                    g = rounds(N)
                    found.append((N, g, s, q))
                    print('found: %d needs %d rounds (%d * %d %s 1)' % (N, g, s, q, '+' if N == s * q + 1 else '-'))
                    sys.stdout.flush()
    if len(FAC) > 200000:
        FAC.clear()
        MEMO.clear()
print('stage 2: %d candidates s q +- 1 with s <= %d, %d primes, %d need 7 rounds or more' % (tried, S2, primes, len(found)))
if found:
    N, g, s, q = min(found)
    print('smallest found: %d (%d rounds)' % (N, g))
print('done')
