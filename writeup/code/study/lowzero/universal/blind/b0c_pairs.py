"""replace the product pairs: the first 10 pairs (SEED+7) were read by accident when bt_test.py was run on the blind folder
before registration (their zero lists exist from round 144): D12 = 11.60+0.84i +- 13.11 (seen, not blind). New pairs: the
138 round-144 keys minus the 20 curves of the seen pairs and the 4 proxy curves, sorted, permuted by default_rng(SEED + 8);
the first 20 form 10 pairs (consecutive)."""
import os, json, numpy as np
SEED = 20261125
d = json.load(open('testset.json'))
z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..', 'data', 'study', 'r144', 'fresh_zeros.npz'))
used = set(k for p in d['T5pairs']['test'] + d['T5pairs']['proxy'] for k in p)
rest = sorted(k for k in z.keys() if k not in used)
rng = np.random.default_rng(SEED + 8); perm = [rest[i] for i in rng.permutation(len(rest))]
new = [[perm[2 * i], perm[2 * i + 1]] for i in range(10)]
d['T5pairs_seen'] = dict(pairs=d['T5pairs']['test'], D12_seen=dict(mean=[11.60, 0.84], se=13.11, decision='PASS'),
                         note='read by accident before registration (bt_test.py run on the blind folder); not blind; reported only')
d['T5pairs']['test'] = new
d['T5pairs']['rule_replacement'] = __doc__
json.dump(d, open('testset.json', 'w'), indent=1)
print('new pairs:', new)
