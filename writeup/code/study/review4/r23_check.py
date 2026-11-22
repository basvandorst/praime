"""Independent recomputation of round 23: N_law(4e18; g) = N(1e6; g) + int_{1e6}^{4e18} sum_j (-1)^j U_j / ln^(j+2) t dt."""
import json, math, sys
import numpy as np
import mpmath as mp
sys.path.insert(0, '/home/bas/scratch/study-review4')
from hl_exact import U_sums

mp.mp.dps = 30
X0, X = 10 ** 6, 4 * 10 ** 18
# counts from OeS table (read raw text)
N = {}
for line in open('/tmp/claude-1000/-home-bas-praime/cfd838fc-8861-4900-931a-268c5006ce39/scratchpad/dl_oes/t0.txt'):
    if line.startswith('#') or not line.strip() or line.startswith('-') or line.strip().startswith('g'):
        continue
    f = line.split()
    g = int(f[0].rstrip('*'))
    N[g] = int(f[2])
# primes to 1e6 + small
lim = 10 ** 6 + 100
s = np.ones(lim + 1, bool); s[:2] = False
for i in range(2, int(lim ** .5) + 1):
    if s[i]:
        s[i * i::i] = False
P = np.nonzero(s)[0]
P = P[: np.searchsorted(P, X0) + 1]       # include first prime above 1e6
assert P[-2] < X0 < P[-1]
d0 = np.diff(P)
I = {}


def Ik(k):
    if k not in I:
        I[k] = mp.quad(lambda u: mp.e ** u / u ** k, mp.linspace(mp.log(X0), mp.log(X), 8))
    return I[k]


rep = json.load(open('/home/bas/praime/writeup/data/study/r23.json'))
out = {}
tot_sets = 0
zs = []
for g in range(2, 72, 2):
    U, cnt = U_sums(g)
    tot_sets += cnt
    law = sum((-1) ** j * U[j] * Ik(j + 2) for j in range(len(U))) + int(np.sum(d0 == g))
    d = N[g]
    z = (d - law) / mp.sqrt(d)
    zs.append(float(z))
    oes = (d - law) / mp.sqrt(2 * d * mp.log(mp.log(d)))
    out[g] = {'law': mp.nstr(law, 20), 'z': float(z), 'z_oes': float(oes),
              'repo_law': rep[str(g)]['law'], 'law_diff_vs_repo_rel': float((law - rep[str(g)]['law']) / law)}
    print(g, d, mp.nstr(law, 19), 'z %+.2f' % float(z), 'oes %.2f' % float(oes),
          'vs repo %+.1e' % out[g]['law_diff_vs_repo_rel'], flush=True)
zs = np.array(zs)
print('total admissible interior sets (incl. empty) for gaps 2..70:', tot_sets)
print('max|z| %.2f, #|z|>1: %d of %d, mean z^2 %.3f, max|rel| %.2e' % (abs(zs).max(), (abs(zs) > 1).sum(), len(zs), (zs ** 2).mean(),
      max(abs(float((N[g] - mp.mpf(out[g]['law'])) / N[g])) for g in out)))
print('max oes ratio', max(abs(v['z_oes']) for v in out.values()))
json.dump(out, open('/home/bas/scratch/study-review4/r23_check.json', 'w'), indent=1)
