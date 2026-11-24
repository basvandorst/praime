"""RETROSPECTIVE: the rounds 117/120 sets whose count falls short of theta/pi (missing zeros), and gamma without them."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *
f = np.load(os.path.join(HERE, 'influence.npz')); c2 = f['c2']; D2 = [int(d) for d in f['D2']]
z = np.load(os.path.join(HERE, 'r120_100_600.npz')); Dz = [int(d) for d in z['D']]; assert Dz == D2
short = [i for i in range(len(D2)) if c2[i, 0] + 1 < -0.5]
print('sets with mean(N - theta/pi) < -0.5 over 100-600:', [(D2[i], round(c2[i, 0] + 1, 3)) for i in short])
keep = np.setdiff1d(np.arange(len(D2)), short); X = design(D2)
print('gamma all %.4f; without them %.4f' % (gamma(z['S'], z['F'], z['Dg'], X), gamma(z['S'][keep], z['F'][keep], z['Dg'][keep], X[keep])))
