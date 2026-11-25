"""time Rankin-Selberg proxies (not in the test set) through the registered zero routine: the drawn proxy and the coprime pair
with the largest N1 N2 <= 1500 that is not a test pair."""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cypari2
from bz4_lib import zeros
pari = cypari2.Pari(); pari.allocatemem(1100 * 10 ** 6)
TS = json.load(open('testset4.json'))
test = set((tuple(x[1]), tuple(y[1])) for x, y in TS['rs'])
cands = []
for i, x in enumerate([]): pass
p0 = TS['rs_proxies'][0]
print('drawn proxy', p0[0][0], 'x', p0[1][0], flush=True)
for T2 in (60,):
    g, info = zeros(pari, ['rs', p0[0][1], p0[1][1], 1], T2)
    print('RS proxy %dx%d T2=%d: n %d status %s drift %.3f r_mean %.3f c_exp %.3f secs %.1f' % (p0[0][0], p0[1][0], T2, info['n'], info['status'], info['drift'], info['r_mean'], info['c_exp'], info['secs']), flush=True)
big = [[37, [0, 0, 1, -1, 0]], [39, [1, 1, 0, -4, -5]]]   # (37, 39) is not a test pair (checked by conductor pair)
for T2 in (60,):
    g, info = zeros(pari, ['rs', big[0][1], big[1][1], 1], T2)
    print('RS large proxy N1 N2 = %d T2=%d: n %d status %s secs %.1f' % (big[0][0] * big[1][0], T2, info['n'], info['status'], info['secs']), flush=True)
