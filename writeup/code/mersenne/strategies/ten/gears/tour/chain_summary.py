"""Collects the results of the chain of wheels (machine 17) into chain.json and chain.txt.
Every number is read from the output file of the script that computed it."""
import json, re, os
from sympy import isprime
from chain_wheel import certify
from chain_verify import verify

out = {}
lines = []


def say(s):
    lines.append(s)
    print(s)


MX = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127}
lit = []
best = {'p': (0, 0), 'm': (0, 0), 'w': (0, 0)}
for t in range(2, 251):
    cert, nodes = certify(t, stop=30)
    if isinstance(cert, list):
        assert verify(t, cert)
        lit.append(t)
    p = [k for k in lit if isprime(k)]
    m = [k for k in p if k in MX]
    w = [k for k in p if isprime(k - 2) or isprime(k + 2)]
    for key, a in (('p', p), ('m', m), ('w', w)):
        if len(a) > best[key][0]:
            best[key] = (len(a), t)
assert all(isprime(k) for k in lit)
out['row250'] = list(best['p'] + best['m'] + best['w'])
say('chain of wheels, 250 wheels: primes %d (tick %d), Mersenne exponents %d (tick %d), twin primes %d (tick %d); lamps of numbers that are not prime: %d' % (tuple(out['row250']) + (sum(1 for k in lit if not isprime(k)),)))

sm = json.load(open('chain_test_small.json'))
out['small'] = sm
for k, v in sm.items():
    assert v['wrong'] == 0 and v['bad'] == 0
    say('Python, %s, every number 2 to %d: %d primes proved, %d primes without a chain, 0 composites given a chain, 0 chains rejected by the checker' % (k, v['top'], v['proved'], v['undecided']))

big = []
pat = re.compile(r'^\s*(\d+) digits, (.*?)\s*: (\d+) primes, proved (\d+), no chain (\d+), timed out \((\d+) s\) (\d+), checker rejects (\d+), median steps (\S+), median time (\S+) s, longest (\S+) s')
pari = {}
for fn in ('chain_test_big.txt', 'chain_test_big2.log'):
    if not os.path.exists(fn):
        continue
    for ln in open(fn):
        m = pat.match(ln)
        if m:
            g = m.groups()
            row = dict(digits=int(g[0]), mode=g[1], primes=int(g[2]), proved=int(g[3]), no_chain=int(g[4]), tlim=int(g[5]), timeout=int(g[6]), rejects=int(g[7]), steps=g[8], median_s=g[9], max_s=g[10])
            assert row['rejects'] == 0
            big.append(row)
        m = re.match(r'^\s*(\d+) digits: PARI isprime \(a proof\) agrees on (\d+) of (\d+)', ln)
        if m:
            assert m.group(2) == m.group(3)
            pari[int(m.group(1))] = int(m.group(2))
out['big'] = big
out['pari'] = pari
for r in big:
    say('Python, %d digits, %s: proved %d of %d (no chain %d, timed out %d), median steps %s, median time %s s' % (r['digits'], r['mode'], r['proved'], r['primes'], r['no_chain'], r['timeout'], r['steps'], r['median_s']))

adv = json.load(open('chain_adversary.json'))
out['adv'] = adv
assert adv['spsp']['chain'] == 0 and adv['mut']['accepted_not_prime'] == 0 and all(adv[k]['chain'] == 0 for k in adv if k.startswith('pq'))
say('composites: %d strong pseudoprimes to base 2 below %d, chains built 0; %d composites p(2p - 1) of 20 to 60 digits that pass the doubling wheel, chains built 0; %d tampered certificates, accepted for a number that is not prime: 0' % (
    adv['spsp']['count'], adv['spsp']['top'], sum(sum(adv[k].values()) for k in adv if k.startswith('pq')), adv['mut']['made']))

js = []
for ln in open('chain_js_check.txt'):
    m = re.match(r'^(.*?): (\d+) numbers, proved prime (\d+), primes without a chain (\d+), composites given a chain (\d+), chains the Python checker rejects (\d+), composites turned away (\d+), steps by kind (\[.*?\])(?:, median time (\d+) ms, longest (\d+) ms)?', ln)
    g = m.groups()
    assert g[4] == '0' and g[5] == '0'
    js.append(dict(label=g[0], numbers=int(g[1]), proved=int(g[2]), no_chain=int(g[3]), median_ms=g[8], max_ms=g[9]))
    d = re.search(r'(\d+) digits', g[0])
    if d:
        js[-1]['digits'] = int(d.group(1))
out['js'] = js
jp = [r for r in js if 'digits' in r]
out['js_primes'] = sum(r['numbers'] for r in jp)
out['js_proved'] = sum(r['proved'] for r in jp)
out['js_maxdigits'] = max(r['digits'] for r in jp)
say('page code under node: %d random primes of %d to %d digits, proved %d, every certificate accepted by the Python checker' % (out['js_primes'], min(r['digits'] for r in jp), out['js_maxdigits'], out['js_proved']))
pc = open('page_check_17.txt').read()
assert 'machine 17: Lamps lit so far: 669. Mistakes so far: 0. script errors 0' in pc
cp = json.load(open('classpoly.json'))
out['families'] = len(cp)
assert '%d checked, 0 differ' % len(cp) in open('classpoly_check.txt').read()
json.dump(out, open('chain.json', 'w'))
open('chain.txt', 'w').write('\n'.join(lines) + '\n')
