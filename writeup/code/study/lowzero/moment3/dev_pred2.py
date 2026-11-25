import os, sys, math, time, itertools, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tmlib as T
DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
t0 = time.time()
rows = []
for tri in itertools.combinations(DEV, 3):
    F = T.Family(tri)
    pr = F.predictions()
    rows.append({'Ds': tri, 'pred': pr, 'h0': float(F.full.hm[0]), 'h0_no3': float(F.hm_no3[0]),
                 'hm': [float(x) for x in F.full.hm], 'g': [float(x) for x in F.full.g]})
    print(tri, ' '.join('%s %.1f' % (k, v) for k, v in pr.items()), ' A3/prodAij %.4f' % (F.full.hm[0] / F.hm_no3[0]), flush=True)
print('time %.1f s' % (time.time() - t0))
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_pred_triples.json'), 'w'), indent=0)
