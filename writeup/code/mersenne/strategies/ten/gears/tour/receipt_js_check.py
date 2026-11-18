"""Receipts built by the page's code (receipt.js under node) are checked here by receipt.py, and receipts built by receipt.py
(honest and tampered) are checked by the page's code. The two programs share only the description of the method."""
import json, subprocess, sys, os
from gmpy2 import mpz
from receipt import Two, run, build, check, lucas_lehmer

S = sys.argv[1]
LAB = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
node = os.path.join(LAB, 'receipt_node.js') if os.path.exists(os.path.join(LAB, 'receipt_node.js')) else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'receipt_node.js')
ps = '3,5,7,11,13,17,19,23,31,61,89,97,107,127,131,521,523,607,1279,2203,2281,3217,4253,4423,4441,9689'
subprocess.run(['node', node, 'make', ps, S + '/rc_js.json'], check=True)
d = json.load(open(S + '/rc_js.json'))
lines = []
lines.append('SHA-256 written for the page against node\'s own on 6 strings: %d equal' % d['sha'])
acc = same = ll = 0
for r in d['rows']:
    G = Two(r['p'])
    y = (mpz(r['y'][0]), mpz(r['y'][1]))
    mus = [(mpz(a), mpz(b)) for a, b in r['mus']]
    yp, saved, _ = run(G, r['k'])
    same += 1 if yp == y and build(G, r['k'], yp, saved) == mus else 0
    acc += 1 if check(G, y, mus) and r['ok'] else 0
    ll += 1 if r['prime'] == r['ll'] == lucas_lehmer(r['p']) else 0
lines.append('receipts built by the page code for %d exponents (3 to 9689): accepted by the Python checker %d; identical to the receipt Python builds %d; verdict equal to Lucas-Lehmer %d' % (len(d['rows']), acc, same, ll))
inp = []
for p in (61, 127, 521, 1279, 1283, 4423):
    G = Two(p)
    k = 4
    y, saved, _ = run(G, k)
    mus = build(G, k, y, saved)
    ser = lambda v: [str(v[0]), str(v[1])]
    inp.append(dict(p=p, y=ser(y), mus=[ser(m) for m in mus], want=True))
    inp.append(dict(p=p, y=ser(((y[0] + 1) % G.N, y[1])), mus=[ser(m) for m in mus], want=False))
    m2 = list(mus); m2[1] = ((m2[1][0] + 1) % G.N, m2[1][1])
    inp.append(dict(p=p, y=ser(y), mus=[ser(m) for m in m2], want=False))
json.dump(inp, open(S + '/rc_py.json', 'w'))
subprocess.run(['node', node, 'check', S + '/rc_py.json', S + '/rc_py_out.json'], check=True)
res = json.load(open(S + '/rc_py_out.json'))
lines.append('receipts built by Python and checked by the page code: %d honest, accepted %d; %d tampered, accepted %d' % (
    sum(1 for r in res if r['want']), sum(1 for r in res if r['want'] and r['ok']), sum(1 for r in res if not r['want']), sum(1 for r in res if not r['want'] and r['ok'])))
t = [r for r in d['rows'] if r['p'] == 9689][0]
lines.append('page code under node, p = 9689: test %d ms, receipt %d ms, check %s ms (%d turns instead of %d)' % (t['tt'], t['tb'], t['tc'], t['turns'], 9687))
print('\n'.join(lines))
open('receipt_js_check.txt', 'w').write('\n'.join(lines) + '\n')
