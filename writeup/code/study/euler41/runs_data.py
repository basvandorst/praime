"""Round 103 data: windows of r consecutive prime values of f(n) = n^2 + n + 41 for n in [1e3, 1e7), and the longest
run. usage: runs_data.py"""
import json, os, time
import numpy as np
import gmpy2

N0, N1 = 10 ** 3, 10 ** 7
t0 = time.time()
isp = np.fromiter((gmpy2.is_prime(n * n + n + 41) for n in range(N0, N1 + 8)), dtype=bool, count=N1 + 8 - N0)
law = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs_law.json')))
out = {'primes_count': int(isp[:N1 - N0].sum()), 'rows': {}}
win = isp.copy()
for r in range(1, 9):
    if r > 1:
        win = win[:-1] & isp[r - 1:r - 1 + len(win) - 0][:len(win) - 0] if False else None
    w = np.ones(N1 - N0, bool)
    for i in range(r):
        w &= isp[i:i + N1 - N0]
    obs = int(w.sum()); exp = law[str(r)]['expected']
    tol = 3 * exp ** 0.5 + 0.01 * exp
    out['rows'][str(r)] = {'observed': obs, 'expected': exp, 'ratio': obs / exp, 'within_rule': abs(obs - exp) <= tol}
    print(r, 'observed', obs, 'expected %.1f' % exp, 'ratio %.4f' % (obs / exp), 'pass' if abs(obs - exp) <= tol else 'FAIL', flush=True)
# longest run of consecutive prime values in [N0, N1)
x = isp[:N1 - N0].astype(np.int8); best = 0; cur = 0; at = None
d = np.diff(np.concatenate([[0], x, [0]]))
starts = np.nonzero(d == 1)[0]; ends = np.nonzero(d == -1)[0]
L = ends - starts; i = int(np.argmax(L))
out['longest_run'] = {'length': int(L[i]), 'start_n': int(starts[i] + N0)}
out['runs_by_length'] = {str(k): int(np.sum(L == k)) for k in range(1, int(L.max()) + 1)}
out['seconds'] = time.time() - t0
print('longest run', out['longest_run'], 'runs >= 8:', int(np.sum(L >= 8)), round(time.time() - t0), 's')
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs_data.json'), 'w'), indent=1)
