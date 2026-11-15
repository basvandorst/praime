"""Part 5 summary: page transformations with the 14 exponents < 1000 as sources vs all 168 primes < 1000.
Hit = destination (inside the circle 1..1000, y != p) is a Mersenne exponent. Fisher exact test per rule
(one-sided: exponent sources hit more often than prime sources), Bonferroni over the rules."""
import json
from scipy.stats import fisher_exact
d = json.load(open('page142.json'))
rows = d['rows']; out = []
pr = lambda s='': (out.append(s), print(s))
pr(f"{len(rows)} transformations; sources: {d['nSrc']} exponents < 1000 and {d['nPr']} primes < 1000; k = 2, circle 1..1000, outside skipped")
pr(f"pooled: exponent sources {d['totEh']} of {d['totE']} destinations are exponents ({d['totEh'] / d['totE']:.1%}); "
   f"prime sources {d['totPh']} of {d['totP']} ({d['totPh'] / d['totP']:.1%})")
res = []
for r in rows:
    if r['dE'] == 0: continue
    a, b = r['hE'], r['dE'] - r['hE']
    c, e = r['hP'] - r['hE'], (r['dP'] - r['dE']) - (r['hP'] - r['hE'])  # primes that are not exponents
    p = fisher_exact([[a, b], [max(c, 0), max(e, 0)]], alternative='greater')[1] if (a + b) and (c + e) else 1.0
    res.append((p, r, a / (a + b), (c / (c + e)) if (c + e) else float('nan')))
res.sort(key=lambda t: t[0])
nrules = len(res)
pr(f'rules with destinations from exponent sources: {nrules}; Bonferroni threshold 0.05/{nrules} = {0.05 / nrules:.2e}')
pr('top 12 by Fisher p (exponent sources vs the other primes):')
for p, r, fe, fp in res[:12]:
    pr(f"  {r['nr']:3d} {r['name'][:38]:38s} {r['hE']:2d}/{r['dE']:2d} = {fe:5.0%} vs other primes {fp:5.1%}  p = {p:.2e}  hits: {r['hits'][:90]}")
full = [r for p, r, fe, fp in res if r['dE'] and r['hE'] == r['dE']]
pr('rules where every destination of an exponent source is an exponent: ' + ', '.join(f"{r['nr']} {r['name']} ({r['hE']}/{r['dE']})" for r in full))
open('page142.txt', 'w').write('\n'.join(out) + '\n')
