"""Champion at every x (not only powers of ten): the list of x where the set of most common gaps changes, for x up to X,
and a check that 6 is the sole champion for every x in [947, X].  usage: champ_all.py X out.txt"""
import sys
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
s = np.ones(X + 1, bool); s[:2] = False
for i in range(2, int(X ** .5) + 1):
    if s[i]: s[i * i::i] = False
P = np.nonzero(s)[0]
cnt = {}; best = set(); changes = []; sole6_from = None; ok = True
for k in range(1, len(P)):
    g = int(P[k] - P[k - 1]); cnt[g] = cnt.get(g, 0) + 1
    m = max(cnt.values()); b = {d for d, v in cnt.items() if v == m}
    if b != best:
        best = b; changes.append((int(P[k]), sorted(b)))
    if P[k] >= 947 and b != {6}: ok = False
s1 = 'changes of the champion set up to %d (x = the prime where it changes): %s' % (X, changes)
s2 = '6 is the sole champion at every x in [947, %d]: %s' % (X, ok)
for t in (s1, s2): print(t); out.write(t + '\n')
