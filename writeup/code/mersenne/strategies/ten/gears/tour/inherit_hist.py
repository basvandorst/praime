"""Third count of the rounds of machine 23: the rule from the top down with PARI (inherit_big.py) for every prime
up to TOP, compared with the histogram of the C program (inherit X), and for the primes of 5 and 6 rounds that the
C run up to 4.2e9 printed.   usage: inherit_hist.py TOP"""
import sys, re, subprocess
import cypari2
from inherit_big import rounds, FAC, MEMO, pari

TOP = int(sys.argv[1])
hist = {}
p = 3
n = 0
while p <= TOP:
    r = rounds(p)
    hist[r] = hist.get(r, 0) + 1
    n += 1
    if len(FAC) > 300000:
        FAC.clear()
        MEMO.clear()
    p = int(pari.nextprime(p + 1))
c = subprocess.run(['./inherit', str(TOP)], capture_output=True, text=True).stdout
ch = {int(a): int(b) for a, b in re.findall(r'depth (\d+): (\d+) primes', c) if int(a) > 0}
print('every prime from 3 to %d, rounds from the top down with PARI: %s' % (TOP, sorted(hist.items())))
print('the C program on the same range:                             %s' % sorted(ch.items()))
print('equal: %s (%d primes)' % (hist == ch, n))
t = open('inherit_c_4e9.out').read()
listed = [(int(a), int(b)) for a, b in re.findall(r'depth (\d) prime (\d+)', t)]
bad = sum(1 for d, q in listed if rounds(q) != d)
print('the %d primes of 5 and 6 rounds printed by the C run up to 4.2e9 (the first 12 of 5 rounds, all 3 of 6 rounds: %s): rounds from the top down differ for %d' % (len(listed), ', '.join(str(q) for d, q in listed if d == 6), bad))
