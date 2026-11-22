"""Round 50 (preregister/round50.txt): the corrected predictability law (r45b) and the information law (r46) in a fresh
window [1e15, 1e15 + 1e7). usage: r50.py"""
import os, json, math, random
import numpy as np
import theory_r3 as R3
import r45, r45b, r46

if __name__ == '__main__':
    L, W = 10 ** 15, 10 ** 7
    d = r45.run(L, W)
    rng = random.Random(50)
    rows = []
    for row in d['rows']:
        ms, how = r45b.mean_S(row['y'], rng)
        pred = min(1.0, row['predicted'] * ms) if row['u'] > 2.0 else 1.0
        ok = abs(row['measured'] - pred) <= 0.002 + 3 * row['sd']
        rows.append({'y': row['y'], 'u': row['u'], 'measured': row['measured'], 'corrected': pred, 'diff': row['measured'] - pred,
                     'sd': row['sd'], 'pass': bool(ok), 'how': how})
        print('y %8d u %5.2f meas %.4f pred %.4f diff %+.4f %s' % (row['y'], row['u'], row['measured'], pred, row['measured'] - pred, '' if ok else 'MISS'), flush=True)
    info = r46.run(L, W, rows)
    for r in info:
        r['pass'] = bool(r['y'] < 100 or abs(r['H_measured'] - r['H_geometric']) <= 0.01)
    out = {'window': [L, W], 'primes': d['primes'], 'hit': rows, 'info': info,
           'hit_misses': [r['y'] for r in rows if not r['pass']], 'info_misses': [r['y'] for r in info if not r['pass']]}
    print('hit misses', out['hit_misses'], 'info misses', out['info_misses'], flush=True)
    json.dump(out, open(os.path.join(R3.OUT, 'r50.json'), 'w'), indent=1)
