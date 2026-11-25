"""Cross part (decade total minus the nine sub-windows) of the remainder R and of X for pairs h = 2..12 on the three
seen decades, against the Hardy-Littlewood cross part; nominal noise 1.33 v / sqrt(K-1).  -> cross_h.json"""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import json, math
import numpy as np
import decomp
HS = (2, 4, 6, 8, 10, 12)
out = {}
for dec, path in (('1e7-1e8', 'win_1e7.npz'), ('1e8-1e9', 'win_1e8.npz'), ('1e9-1e10', 'win_1e9.npz')):
    tot = decomp.analyse(path, HS); subs = [decomp.analyse(path, HS, sub=s) for s in range(9)]
    for h in HS:
        for Q in (210, 2310, 30030):
            k = 'h%d_%d' % (h, Q); K = tot[k]['K']
            ms = np.array([s[k]['pairs'] for s in subs], float); m = ms.sum()
            r = {}
            for key in ('X', 'R'):
                wm = float(np.sum(np.array([s[k]['measured_u2' + key] for s in subs]) * ms) / m)
                wp = float(np.sum(np.array([s[k]['pred_u2' + key] for s in subs]) * ms) / m)
                cm = tot[k]['measured_u2' + key] - wm; cp = tot[k]['pred_u2' + key] - wp
                sd = 2 * wp * math.sqrt((1 - np.sum(ms ** 2) / m ** 2) / 2 / (K - 1))
                r[key] = {'cross_meas': cm, 'cross_HL': cp, 'sd': sd, 'z': (cm - cp) / sd, 'within_ratio': wm / wp}
            out['%s|%s' % (dec, k)] = r
    print(dec, 'done', flush=True)
json.dump(out, open('cross_h.json', 'w'), indent=1)
for Q in (210, 2310, 30030):
    print('Q', Q)
    for h in HS:
        zs = [out['%s|h%d_%d' % (d, h, Q)]['R']['z'] for d in ('1e7-1e8', '1e8-1e9', '1e9-1e10')]
        rat = [out['%s|h%d_%d' % (d, h, Q)]['R']['cross_meas'] / out['%s|h%d_%d' % (d, h, Q)]['R']['cross_HL'] for d in ('1e7-1e8', '1e8-1e9', '1e9-1e10')]
        print('  h=%2d  R cross z: %s  joint %.2f   meas/HL: %s' % (h, ' '.join('%+.2f' % z for z in zs), sum(zs) / math.sqrt(3), ' '.join('%.2f' % x for x in rat)))
