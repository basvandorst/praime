"""Other pairs (n, n+h), h = 4..12, on the three seen decades: measured / predicted spread (X), remainder (R), beta."""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import json
import decomp
out = {}
for name, path in (('1e7-1e8', 'win_1e7.npz'), ('1e8-1e9', 'win_1e8.npz'), ('1e9-1e10', 'win_1e9.npz')):
    out[name] = decomp.analyse(path, (4, 6, 8, 10, 12))
    for h in (4, 6, 8, 10, 12):
        print(name, h, '  '.join('Q%d X %.2f R %.2f b %.2f' % (Q, out[name]['h%d_%d' % (h, Q)]['ratio_X'], out[name]['h%d_%d' % (h, Q)]['ratio_R'],
              out[name]['h%d_%d' % (h, Q)]['beta_hat']) for Q in (210, 2310, 30030)), flush=True)
json.dump(out, open('h_decomp.json', 'w'), indent=1)
