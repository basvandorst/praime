import numpy as np, math
from q import Qparts, Z
_, s, o, _ = Qparts(Z); Qall = 2*s + 2*o
for M in (3, 5, 10, 30, 100, 300, 1000):
    _, s2, o2, _ = Qparts(Z[:M]); print('first', M, 'zeros (gamma<=%.0f): fraction of Q' % Z[M-1], np.round((2*s2+2*o2)/Qall, 3))
