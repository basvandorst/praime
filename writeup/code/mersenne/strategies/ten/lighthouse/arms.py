"""Prime lighthouse: in the polar plot (r, theta) = (p, p radians) the primes line up on arms because 2 pi ~ 44/7 and
~ 710/113; the arm of p is p mod 44 (or mod 710).  Arms with residue sharing a factor with q hold at most one prime, so
there are phi(44) = 20 and phi(710) = 280 arms.  Counts per arm against Dirichlet (equal), chi-square, and the race
between square and non-square residues (Chebyshev bias: squares mod q lag).  Prediction written before the run: the
arms are equal within chi-square, the square residues are behind on average.   usage: arms.py X out.txt"""
import sys, math
import numpy as np
from scipy.stats import chi2
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
pr('X = %d; 2 pi = %.6f, 44/7 = %.6f, 710/113 = %.6f' % (X, 2 * math.pi, 44 / 7, 710 / 113))
for q in (44, 710):
    c = np.bincount(P % q, minlength=q); arms = [r for r in range(q) if math.gcd(r, q) == 1]; a = c[arms]
    big = int(np.sum(c > 1)); e = a.mean(); x2 = float(np.sum((a - e) ** 2 / e)); dof = len(arms) - 1
    sq = set((t * t) % q for t in range(q) if math.gcd(t, q) == 1)
    s_sq = a[[r in sq for r in arms]].mean(); s_ns = a[[r not in sq for r in arms]].mean()
    pr('q = %d: arms with more than one prime %d (phi = %d); counts min %d max %d mean %.1f; chi2 = %.1f on %d dof, p = %.3f' % (q, big, len(arms), a.min(), a.max(), e, x2, dof, chi2.sf(x2, dof)))
    pr('   %d square residues: mean count %.1f; %d non-squares: %.1f; difference %.1f' % (len(sq), s_sq, len(arms) - len(sq), s_ns, s_ns - s_sq))
