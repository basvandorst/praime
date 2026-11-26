import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P
from collections import Counter
H = os.path.dirname(os.path.abspath(__file__))
for r2 in (0, 1, 2):
    z = np.load(os.path.join(H, 'data', 'disc_r%d.npz' % r2)); ty = z['types']; n = len(ty)
    print('r2=%d (%d fields): type  freq_obs (count)  mass/Z  [obs/mass]' % (r2, n))
    for k, l in enumerate((2, 3, 5, 7)):
        c = Counter(ty[:, k]); Z = 1 + 1 / l + 2 / l**2 + 1 / l**3
        unk = set(c) - set(P.TYPES)
        row = ['%s %.3f/%.3f' % (nm, c.get(nm, 0) / n, P.TYPES[nm][0](l) / Z) for nm in P.TYPES]
        print('  l=%d: ' % l + ', '.join(row) + ('  UNKNOWN %s' % unk if unk else ''))
