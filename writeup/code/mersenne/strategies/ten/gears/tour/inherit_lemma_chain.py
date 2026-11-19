"""The lemma of machine 23 used twice (see inherit_lemma.py): from every prime of 5 or more rounds up to B to
every prime of 6 or more rounds up to 4 B, and from those to every prime of 7 or more rounds up to 16 B.
Why 4 B: a prime q <= 4 B of g + 1 rounds has a helper r of g or more rounds in q - 1 = s r and one in
q + 1 = s' r'. If both were above B, then s = s' = 2 and r' = r + 1, two primes above 2 next to each other.
So q = s r + 1 or s r - 1 with r <= B.
usage: inherit_lemma_chain.py B LIST    LIST holds every prime of 5 or more rounds up to B (lines: rounds prime)"""
import sys
from inherit_big import rounds, atmost, pari

B = int(sys.argv[1])
deep5 = [int(l.split()[1]) for l in open(sys.argv[2]) if int(l.split()[1]) <= B]


def up(helpers, g, top):
    """Every prime up to top of g + 1 or more rounds, from the list of all primes of g or more rounds up to top / 4."""
    out = set()
    tried = 0
    for r in helpers:
        s = 2
        while s * r - 1 <= top:
            for q in (s * r - 1, s * r + 1):
                if q <= top:
                    tried += 1
                    if pari.ispseudoprime(q) and not atmost(q, g):
                        out.add(q)
            s += 2
    return sorted(out), tried


six, t6 = up(deep5, 5, 4 * B)
print('from the %d primes of 5 or more rounds up to %d: %d candidates s r +- 1 up to %d, primes of 6 or more rounds: %d (most rounds %d)' % (len(deep5), B, t6, 4 * B, len(six), max(rounds(q) for q in six)))
seven, t7 = up(six, 6, 16 * B)
print('from these: %d candidates s q +- 1 up to %d, primes of 7 or more rounds: %s' % (t7, 16 * B, ', '.join(str(x) for x in seven) or 'none'))
