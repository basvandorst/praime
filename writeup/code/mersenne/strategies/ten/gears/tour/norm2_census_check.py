"""The census of starts of norm 2 over the whole list (norm2.c census 8, two halves): totals per start,
and every sampled line recomputed in plain Python (Jacobi symbol and both cells after n turns)."""
import json, re
from norm2 import jacobi, power

T = json.load(open('tournament.json'))
starts = {}
lines = last = mine = 0
for i in (0, 1):
    e = open('norm2_census_%d.err' % i).read()
    assert 'done' in e
    g = re.search(r'lines read (\d+), last (\d+), mine (\d+)', e).groups()
    lines, last = int(g[0]), int(g[1])
    mine += int(g[2])
    for a, c, mt, mh, ft, fp, fm, sh in re.findall(r'start (\d+) \+ sqrt (\d+): does not fit (\d+), of these home (\d+) \| fits (\d+), of these \(a \+ sqrt c\)\^n = a \+ sqrt c (\d+), = a - sqrt c (\d+) \| shares a factor (\d+)', e):
        d = starts.setdefault((int(a), int(c)), dict(a=int(a), c=int(c), misfit=0, misfit_home=0, fit=0, fit_same=0, fit_swapped=0, shares=0))
        d['misfit'] += int(mt); d['misfit_home'] += int(mh); d['fit'] += int(ft); d['fit_same'] += int(fp); d['fit_swapped'] += int(fm); d['shares'] += int(sh)
assert 'pipeline exit 0' in open('norm2_census_1.err').read()
sha = open('norm2_census.sha').read().split()[0]
assert sha == '3b92a884be4a365ca91cbd7c9272620bc7e081335a0251e3531b404f2d36f8ef'
assert lines == T['list_count'] == mine
samples = bad = 0
homes = []
for i in (0, 1):
    for l in open('norm2_census_%d.out' % i):
        p = l.split()
        if p[0] == 'misfit':
            homes.append(l.strip())
        elif p[0] == 'sample':
            n, a, c, j, x, y = [int(v) for v in p[1:]]
            samples += 1
            if jacobi(c, n) != j or power(a, c, n, n) != (x, y):
                bad += 1
rows = [starts[k] for k in sorted(starts)]
for d in rows:
    assert d['misfit'] + d['fit'] + d['shares'] == lines
tot = dict(misfit=sum(d['misfit'] for d in rows), misfit_home=sum(d['misfit_home'] for d in rows), fit=sum(d['fit'] for d in rows), fit_same=sum(d['fit_same'] for d in rows), fit_swapped=sum(d['fit_swapped'] for d in rows))
txt = ['census over the list: %d numbers on %d starts; sha256 %s' % (lines, len(rows), sha)]
for d in rows:
    txt.append('start %d + sqrt %d: does not fit %d, home %d | fits %d, cells unchanged %d, cells swapped %d | shares a factor %d' % (d['a'], d['c'], d['misfit'], d['misfit_home'], d['fit'], d['fit_same'], d['fit_swapped'], d['shares']))
txt.append('all starts: tests on a start the number does not fit %d, came home %d; tests on a start it fits %d, cells unchanged %d, cells swapped %d' % (tot['misfit'], tot['misfit_home'], tot['fit'], tot['fit_same'], tot['fit_swapped']))
txt.append('sampled results recomputed in plain Python: %d, different: %d' % (samples, bad))
assert bad == 0 and len(homes) == tot['misfit_home']
print('\n'.join(txt))
open('norm2_census.txt', 'w').write('\n'.join(txt) + '\n')
json.dump(dict(lines=lines, rows=rows, total=tot, samples=samples, sample_bad=bad, homes=homes), open('norm2_census.json', 'w'))
