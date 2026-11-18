"""Collects the results of the two-cell wheel from a start of norm 2 (machine 19) into norm2.json and norm2.txt.
Every number is read from the output file of the program that computed it."""
import json, re
from sympy import isprime
from norm2 import home

out = {}
lines = []


def say(s):
    lines.append(s)
    print(s)


MX = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127}
lit = []
best = {'p': (0, 0), 'm': (0, 0), 'w': (0, 0)}
for t in range(2, 251):
    if home(t):
        lit.append(t)
    p = [k for k in lit if isprime(k)]
    m = [k for k in p if k in MX]
    w = [k for k in p if isprime(k - 2) or isprime(k + 2)]
    for key, a in (('p', p), ('m', m), ('w', w)):
        if len(a) > best[key][0]:
            best[key] = (len(a), t)
assert all(isprime(k) for k in lit)
out['row250'] = list(best['p'] + best['m'] + best['w'])
say('norm-2 wheel, 250 wheels: primes %d (tick %d), Mersenne exponents %d (tick %d), twin primes %d (tick %d); lamps of numbers that are not prime: %d' % (tuple(out['row250']) + (sum(1 for k in lit if not isprime(k)),)))
L = json.load(open('norm2_list.json'))
assert L['home'] == 0 and L['sample_bad'] == 0 and L['lines'] == L['list_count']
out['list'] = L
say('list below 2^64: %d numbers, %d on a wheel, %d came home, %d turned away for a common factor with a small step, squares %s, largest a %d; %d sampled lines recomputed in Python, %d differ' % (L['lines'], L['tested'], L['home'], L['factor'], L['squares'], L['amax'], L['samples'], L['sample_bad']))
C = open('norm2_check.txt').read()
m = re.search(r'base-2 pseudoprimes below 10\^10: (\d+); came home on the norm-2 wheel: (\d+)', C).groups()
m2 = re.search(r'odd composites below (\d+) on wheels they fit .*?: (\d+) cases, first .*?; of these not a base-2 pseudoprime: (\d+)', C).groups()
assert m[1] == '0' and m2[2] == '0'
out['lemma'] = dict(top=int(m2[0]), cases=int(m2[1]), psp1e10=int(m[0]))
say('norm lemma where it is not empty: %d cases below %d, all base-2 pseudoprimes' % (int(m2[1]), int(m2[0])))
rng = []
for part in 'ab':
    e = open('norm2_range_%s.err' % part).read()
    g = re.search(r'norm-2 wheel, odd n in \[(\d+), (\d+)\): came home (\d+), perfect squares (\d+), sharing a factor with a small c (\d+), largest a among those home (\d+) \(n = (\d+)\)', e).groups()
    assert e.strip().endswith('done')
    rng.append([int(x) for x in g])
P = open('norm2_pcount.txt').read()
pc = {(int(a), int(b)): int(c) for a, b, c in re.findall(r'odd primes in \[(\d+), (\d+)\): (\d+)', P)}
tot_home = tot_pr = 0
for lo, hi, hm, sq, fac, amax, amaxn in rng:
    assert pc[(lo, hi)] == hm, (lo, hi)
    tot_home += hm
    tot_pr += pc[(lo, hi)]
out['range'] = dict(lo=rng[0][0], hi=rng[-1][1], home=tot_home, primes=tot_pr, amax=max(r[5] for r in rng), amaxn=[r[6] for r in rng if r[5] == max(q[5] for q in rng)][0])
say('every odd number in [%d, %d), no list: came home %d; odd primes there by an independent sieve %d; largest a for a prime %d (n = %d)' % (out['range']['lo'], out['range']['hi'], tot_home, tot_pr, out['range']['amax'], out['range']['amaxn']))
B = open('norm2_bench.txt').read()
bench = {}
for part, key in (('primes', 'primes'), ('odd numbers', 'odd')):
    blk = B.split('== %s' % part)[1].split('==')[0]
    d = {}
    for name, cnt, yes, ns in re.findall(r'^(.*?): (\d+) numbers, said prime (\d+), (\d+) ns per number', blk, re.M):
        d[name.split(' (')[0]] = dict(count=int(cnt), yes=int(yes), ns=int(ns))
    assert len(set(v['yes'] for v in d.values())) == 1
    bench[key] = d
out['bench'] = bench
for key in bench:
    say('speed on 64-bit %s: ' % key + '; '.join('%s %d ns' % (k, v['ns']) for k, v in bench[key].items()))
assert 'machine 19: Lamps lit so far: 669. Mistakes so far: 0. script errors 0' in open('page_check_19.txt').read()
json.dump(out, open('norm2.json', 'w'))
open('norm2.txt', 'w').write('\n'.join(lines) + '\n')
