# Puzzle 1281 (writeup/1281.md): the largest prime whose blocks of k adjacent digits are pairwise coprime k-digit
# numbers (case 1: k = 2; case 3: k = 2 and no block is prime; case 4: k = 3 and no block is prime).
# Digits are 1..9: a block starting with 0 is not a k-digit number, and a 0 anywhere starts a block (except as the
# last digit, where the number would not be prime).
# Exhaustive depth-first search over digit strings: state = (last k-1 digits, set of primes used by the blocks so
# far). Blocks are pairwise coprime iff their sets of prime factors are disjoint, so a block may be appended iff none
# of its primes is in use. For every length the search keeps the largest prime string (primality checked only when
# the string beats the current best of its length), so the output gives the largest prime of every possible length.
# Usage: python3 p1281_search.py k noprimes(0/1)
import sys
from sympy import factorint, isprime

k = int(sys.argv[1]); noprimes = sys.argv[2] == '1'
blocks = {}
for b in range(10 ** (k - 1), 10 ** k):
    s = str(b)
    if '0' in s: continue
    f = factorint(b)
    if noprimes and sum(f.values()) == 1: continue
    blocks[s] = f
primes = sorted({p for f in blocks.values() for p in f}); pid = {p: i for i, p in enumerate(primes)}
mask = {s: sum(1 << pid[p] for p in f) for s, f in blocks.items()}
best = {}                      # length -> largest prime string of that length
nodes = [0]; longest = [0]
def dfs(s, used):
    nodes[0] += 1
    L = len(s); longest[0] = max(longest[0], L)
    if s[-1] in '1379' and (L not in best or s > best[L]) and isprime(int(s)): best[L] = s
    for d in '123456789':
        t = s[-(k - 1):] + d
        m = mask.get(t)
        if m is not None and not (m & used): dfs(s + d, used | m)
for s, m in mask.items(): dfs(s, m)
print('k=%d%s: %d blocks over %d primes, %d strings visited, longest string %d digits' %
      (k, ', no prime blocks' if noprimes else '', len(blocks), len(primes), nodes[0], longest[0]))
for L in sorted(best, reverse=True)[:4]: print('  largest prime with %d digits: %s' % (L, best[L]))
print('answer (largest prime): %s' % best[max(best)])
