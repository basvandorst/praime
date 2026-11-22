"""Round 112 check: the true regions up to 1e9 where pi(x;4,1) > pi(x;4,3) (segmented sieve), the share of primes there,
and the zero prediction's share restricted to x <= 1e9. usage: r112_check.py"""
import os, sys, json, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'study'))
import common as C
H = os.path.dirname(os.path.abspath(__file__))
X = 10 ** 9; SEG = 10 ** 7; lead = 0; regions = []; inside = False; start = None; nprimes = 0; nlead = 0
for s0 in range(0, X, SEG):
    n = min(SEG, X - s0)
    P = np.nonzero(C.primes_mask(max(s0, 3), n - (3 if s0 == 0 else 0)))[0].astype(np.int64) + max(s0, 3)
    step = np.where(P % 4 == 1, 1, -1)
    run = lead + np.cumsum(step)                     # pi(p;4,1) - pi(p;4,3) after each prime p
    pos = run > 0
    nprimes += len(P); nlead += int(pos.sum())
    ch = np.nonzero(np.diff(np.concatenate([[inside], pos]).astype(int)))[0]
    for c in ch:
        if not inside: start = int(P[c]); inside = True
        else: regions.append((start, int(P[c]))); inside = False
    lead = int(run[-1]) if len(run) else lead
if inside: regions.append((start, X))
d = json.load(open(os.path.join(H, 'r112.json')))
out = {'X': X, 'share_primes_class1_ahead': nlead / nprimes, 'regions': regions[:20], 'n_regions': len(regions)}
print('true share to 1e9: %.5f; %d regions; first: %s' % (nlead / nprimes, len(regions), regions[:8]))
json.dump(out, open(os.path.join(H, 'r112_check.json'), 'w'), indent=1)
