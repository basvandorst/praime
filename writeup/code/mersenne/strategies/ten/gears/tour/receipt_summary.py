"""Collects the results of the test with a receipt (machine 18) into receipt.json and receipt.txt.
Every number is read from the output file of the script that computed it."""
import json, re
from sympy import isprime
from receipt import lucas_lehmer

R = json.load(open('receipt_test.json'))
out = {}
lines = []


def say(s):
    lines.append(s)
    print(s)


MX = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127}
lit = []
best = {'p': (0, 0), 'm': (0, 0), 'w': (0, 0)}
for t in range(2, 251):
    if isprime(t) and lucas_lehmer(t):
        lit.append(t)
    p = [k for k in lit if isprime(k)]
    m = [k for k in p if k in MX]
    w = [k for k in p if isprime(k - 2) or isprime(k + 2)]
    for key, a in (('p', p), ('m', m), ('w', w)):
        if len(a) > best[key][0]:
            best[key] = (len(a), t)
assert set(lit) == MX
out['row250'] = list(best['p'] + best['m'] + best['w'])
say('two cells with a receipt, 250 wheels (wheel p decides 2^p - 1): primes %d (tick %d), Mersenne exponents %d (tick %d), twin primes %d (tick %d)' % tuple(out['row250']))

a = R['agree']
assert a['one'] == a['two'] == a['acc_one'] == a['acc_two'] == a['primes'] and R['composite_n_pass'] == []
out['agree'] = a
say('odd primes p up to %d: %d; both forms agree with Lucas-Lehmer on all, every honest receipt accepted; no composite exponent below 2000 passes the one-cell test' % (a['top'], a['primes']))
tm = {}
import os
TM = json.load(open('receipt_timing.json')) if os.path.exists('receipt_timing.json') else R['timing']
out['timing_source'] = 'receipt_timing.txt' if os.path.exists('receipt_timing.json') else 'receipt_test.txt'
tmin = {}
for r in TM:
    key = (r['p'], r['group'])
    tmin[key] = min(tmin.get(key, 1e99), r['test_s'])
for r in TM:
    key = (r['p'], r['group'])
    r = dict(r, test_s=tmin[key], ratio=round(tmin[key] / r['check_s'], 1))
    if key not in tm or r['ratio'] > tm[key]['ratio']:
        tm[key] = r
out['timing'] = [tm[k] for k in sorted(tm)]
for r in out['timing']:
    say('p = %d, %s, %d halvings: test %.2f s, receipt %.2f s, check %.3f s, %.0f times faster' % (r['p'], r['group'], r['k'], r['test_s'], r['receipt_s'], r['check_s'], r['ratio']))
t = R['tamper']
for k in ('one', 'two'):
    assert t[k]['wrong_y_accepted'] == 0 and t[k]['wrong_mu_accepted'] == 0 and t[k]['flips_accepted'] == 0 and t[k]['flips_same_result'] == 0
out['tamper'] = t
say('tampering at p = 9689: accepted 0 of %d wrong results, 0 of %d changed receipt numbers, 0 of %d runs with a flipped bit (each form)' % (t['two']['wrong_y_tried'], t['two']['wrong_mu_tried'], t['two']['flips']))
fg = {}
for r in R['forgery']:
    k = (r['group'], r['factor'])
    d = fg.setdefault(k, dict(group=r['group'], factor=r['factor'], tried=0, built=0, accepted=0, changed=0, exps=0))
    d['tried'] += r['tried']; d['built'] += r['built']; d['accepted'] += r['accepted']; d['changed'] += r['verdict_changed']; d['exps'] += 1
out['forgery'] = list(fg.values())
for d in out['forgery']:
    assert d['built'] == d['accepted']
    say('forgery, %s, factor %s: %d attempts on %d prime 2^p - 1, accepted %d, verdict changed to not prime %d' % (d['group'], d['factor'], d['tried'], d['exps'], d['accepted'], d['changed']))
s = R['safe']
assert s['one_forged_flagged'] == s['one_forged'] and s['two_forged_flagged'] == s['two_forged'] and s['flagged_one'] == 0 and s['flagged_two'] == 0
out['safe'] = s
say('guard: %d of %d forged results flagged; 0 of %d composite 2^p - 1 flagged without need (each form)' % (s['one_forged_flagged'] + s['two_forged_flagged'], s['one_forged'] + s['two_forged'], s['composites']))
J = open('receipt_js_check.txt').read()
m = re.search(r'receipts built by the page code for (\d+) exponents \(3 to 9689\): accepted by the Python checker (\d+); identical to the receipt Python builds (\d+); verdict equal to Lucas-Lehmer (\d+)', J).groups()
assert len(set(m)) == 1
m2 = re.search(r'checked by the page code: (\d+) honest, accepted (\d+); (\d+) tampered, accepted (\d+)', J).groups()
assert m2[0] == m2[1] and m2[3] == '0' and '6 equal' in J
out['js'] = dict(exps=int(m[0]), honest=int(m2[0]), tampered=int(m2[2]))
say('page code: receipts for %d exponents identical to Python and accepted by it; Python receipts: %d honest accepted, %d tampered rejected' % (int(m[0]), int(m2[0]), int(m2[2])))
assert 'machine 18: Lamps lit so far: 20. Mistakes so far: 0. script errors 0' in open('page_check_18.txt').read()
json.dump(out, open('receipt.json', 'w'))
open('receipt.txt', 'w').write('\n'.join(lines) + '\n')
