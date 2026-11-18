"""Second check of the run over the list: every sampled line of norm2_list.out is recomputed in plain Python (norm2.py),
and the count of lines, the last number and the checksum of the download are compared with the first run of section 20."""
import json, re
from norm2 import start, power, home

T = json.load(open('tournament.json'))
err = open('norm2_list.err').read()
m = re.search(r'lines (\d+), last (\d+), on a wheel (\d+), came home (\d+), perfect squares (\d+), sharing a factor with a small c (\d+), largest a (\d+) \(n = (\d+)\)', err)
lines, last, tested, came, sq, fac, amax, amaxn = [int(x) for x in m.groups()]
assert 'pipeline exit 0' in err and err.count('done') == 1
sha = open('norm2_list.sha').read().split()[0]
out = open('norm2_list.out').read().split('\n')
squares = [int(l.split()[1]) for l in out if l.startswith('square')]
homes = [l for l in out if l.startswith('home')]
samples = [l.split() for l in out if l.startswith('sample')]
bad = 0
for s in samples:
    n, a, c, x, y, r = [int(v) for v in s[1:]]
    st = start(n)
    if st is None:
        ok = r == -1
    elif st[0] == 'factor':
        ok = r == -2
    else:
        ok = st == (a, c) and power(a, c, n, n) == (x, y) and (r == 1) == home(n)
    bad += 0 if ok else 1
res = dict(lines=lines, last=last, tested=tested, home=came, squares=squares, factor=fac, amax=amax, amaxn=amaxn, sha256=sha, samples=len(samples), sample_bad=bad, list_count=T['list_count'])
txt = ['lines read from the list: %d (the earlier runs of section 20 read %d); last number %d' % (lines, T['list_count'], last),
       'sha256 of the downloaded file: %s' % sha,
       'put on a norm-2 wheel: %d; came home: %d; perfect squares (no wheel, so not prime): %s; turned away because they share a factor with a small step: %d; largest a: %d (n = %d)' % (tested, came, squares, fac, amax, amaxn),
       'every 4,000th line recomputed in plain Python (start, both cells after n turns, verdict): %d lines, different: %d' % (len(samples), bad)]
assert lines == T['list_count'] and tested + len(squares) + fac == lines and came == 0 and not homes and bad == 0
assert sha == '3b92a884be4a365ca91cbd7c9272620bc7e081335a0251e3531b404f2d36f8ef'
import os
if os.path.exists('norm2_second_0.err') and os.path.exists('norm2_second_1.err'):
    tot = dict(mine=0, tested=0, home=0, sq=0, fac=0, amax=0)
    for i in (0, 1):
        e = open('norm2_second_%d.err' % i).read()
        g = re.search(r'lines read (\d+), last (\d+), mine (\d+), on a wheel (\d+), came home (\d+), perfect squares (\d+), sharing a factor with a small c (\d+), largest a (\d+)', e).groups()
        assert 'done' in e and int(g[0]) == lines and int(g[1]) == last and 'no step' not in e
        tot['mine'] += int(g[2]); tot['tested'] += int(g[3]); tot['home'] += int(g[4]); tot['sq'] += int(g[5]); tot['fac'] += int(g[6]); tot['amax'] = max(tot['amax'], int(g[7]))
    assert 'pipeline exit 0' in open('norm2_second_1.err').read()
    assert open('norm2_second.sha').read().split()[0] == sha
    assert tot['mine'] == lines and tot['tested'] == tested and tot['home'] == 0 and tot['sq'] == len(squares) and tot['fac'] == fac and tot['amax'] == amax
    assert open('norm2_second_0.out').read().strip() == '' and open('norm2_second_1.out').read().strip() == ''
    res['second'] = tot
    txt.append('second program (norm2_second.c, written separately: plain remainders, power from the lowest bit up), the whole list again in two halves: %d numbers, on a wheel %d, came home %d, perfect squares %d, turned away for a common factor %d, largest a %d: every count equal to the first run' % (tot['mine'], tot['tested'], tot['home'], tot['sq'], tot['fac'], tot['amax']))
print('\n'.join(txt))
open('norm2_list_check.txt', 'w').write('\n'.join(txt) + '\n')
json.dump(res, open('norm2_list.json', 'w'))
