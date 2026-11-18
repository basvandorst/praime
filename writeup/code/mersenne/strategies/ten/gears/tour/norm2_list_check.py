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
print('\n'.join(txt))
open('norm2_list_check.txt', 'w').write('\n'.join(txt) + '\n')
json.dump(res, open('norm2_list.json', 'w'))
