import json, sys, time
import numpy as np
from common import make_nulls, A, L52

n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
t = time.time()
na, nb = make_nulls(n)
json.dump({'real': A, 'a': na, 'b': nb}, open('nulls.json', 'w'))
for name, S in (('a', na), ('b', nb)):
    S = np.array(S, dtype=float)
    l = np.log2(S)
    print(name, 'n=%d' % len(S), 'log2 a52 mean %.2f sd %.2f' % (l[:, -1].mean(), l[:, -1].std()),
          'a1..5 example', [int(v) for v in S[0, :8]],
          'count<=127 mean %.1f (real 12)' % (S <= 127).sum(1).mean(),
          'count<=1e4 mean %.1f (real %d)' % ((S <= 10 ** 4).sum(1).mean(), sum(x <= 10 ** 4 for x in A)))
print('real log2 a52 %.2f' % L52, 'time %.1fs' % (time.time() - t))
