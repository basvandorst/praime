# Turn a float optimum (s2_opt_*.json) into exact rational measures for the certifier.
import json, sys
from fractions import Fraction as F

d = json.load(open(sys.argv[1]))
out = sys.argv[2]
K = d['K']
DEN_R = 4096          # free radii rounded to multiples of 1/4096
DEN_W = 2**20         # weights rounded to multiples of 2^-20
radii, mu, nu = d['radii'], d['mu'], d['nu']
atoms = {}
for i, (r, m, n) in enumerate(zip(radii, mu, nu)):
    rr = F(i + 1, K) if i < K else F(round(r*DEN_R), DEN_R)
    wm, wn = F(round(m*DEN_W), DEN_W), F(round(n*DEN_W), DEN_W)
    if wm == 0 and wn == 0:
        continue
    pm, pn = atoms.get(rr, (F(0), F(0)))
    atoms[rr] = (pm + wm, pn + wn)
rs = sorted(atoms)
mus = [atoms[r][0] for r in rs]
nus = [atoms[r][1] for r in rs]
# exact normalisation: put the rounding defect on the largest weight
for w in (mus, nus):
    i = max(range(len(w)), key=lambda k: w[k])
    w[i] += 1 - sum(w)
assert sum(mus) == 1 and sum(nus) == 1 and min(mus) >= 0 and min(nus) >= 0 and 0 < rs[0] and rs[-1] <= 1
json.dump({'radii': [str(r) for r in rs], 'mu': [str(m) for m in mus], 'nu': [str(n) for n in nus],
           'source': sys.argv[1]}, open(out, 'w'), indent=1)
print(len(rs), 'atoms', [float(r) for r in rs])
print('mu', [round(float(m), 5) for m in mus])
print('nu', [round(float(n), 5) for n in nus])
