"""Folded wheels alone (no doubling wheel, no list): every odd number below the bound that is not prime and does not fit the ring
was put on the smallest fold that says something and on which n is not in H (fold.c -show). Which folds make mistakes?
A pass is listed with the order of n on the fold: the smallest k with n^k in H."""
import re, sys, json
from fold_check import prim_root, fold_vec
tag = sys.argv[1] if len(sys.argv) > 1 else '1e7'
err = open('folded_alone_%s.err' % tag).read(); out = open('folded_alone_%s.out' % tag).read().split('\n')
assert not any(l.startswith('ring') for l in out)
P = {}
for l in out:
    f = l.split()
    if len(f) == 6: P.setdefault((int(f[1]), int(f[2])), []).append(int(f[5]))
def order(n, r, e):
    g = prim_root(r); ind = {}; x = 1
    for k in range(r-1): ind[x] = k; x = x*g % r
    i = ind[n % r] % e; k = 1
    while (i*k) % e: k += 1
    return k
from fold_check import folds_of
def orders_tested(r, e):
    g = prim_root(r); ind = {}; x = 1
    for k in range(r-1): ind[x] = k; x = x*g % r
    es = [q for q in folds_of(r) if ind[2] % q != 0]; out = set()
    for x in range(2, r-1):
        pick = next((q for q in es if ind[x] % q != 0), None)
        if pick == e: out.add(order(x, r, e))
    return sorted(out)
R = {'bound': tag, 'folds': []}; lines = []
for m in re.finditer(r'ring (\d+): tested (\d+);(.*?) no fold: (\d+) tested', err):
    r = int(m.group(1))
    for e, a, b in re.findall(r'fold of (\d+) cells: (\d+) tested, (\d+) passed', m.group(3)):
        e, a, b = int(e), int(a), int(b); ps = sorted(P.get((r, e), [])); assert len(ps) == b
        byord = {}
        for n in ps: byord[order(n, r, e)] = byord.get(order(n, r, e), 0) + 1
        chk = sum(1 for n in ps[:40] if fold_vec(n, r, e)[0])
        R['folds'].append({'r': r, 'e': e, 'orders_tested': orders_tested(r, e), 'tested': a, 'passed': b, 'first': ps[:3], 'by_order': byord, 'python_confirms_first': [chk, min(40, len(ps))]})
        lines.append('ring of %d cells folded onto %d (order of n on the fold among the tested numbers: %s): %s tested, %s came home%s%s' % (r, e, ' or '.join(map(str, orders_tested(r, e))), format(a, ','), format(b, ','), (', first ' + ', '.join(map(str, ps[:3]))) if ps else '', ('; order of n on the fold among them: ' + ', '.join('%d: %d' % kv for kv in sorted(byord.items())) + '; plain Python confirms %d of the first %d' % (chk, min(40, len(ps)))) if ps else ''))
    for e in re.findall(r'fold of (\d+) cells: says nothing', m.group(3)): lines.append('ring of %d cells folded onto %s: says nothing (2 is in H)' % (r, e))
hi = [x for x in R['folds'] if any(k >= 3 for k in x['by_order'])]
R['passes_with_order_3_or_more'] = sum(v for x in R['folds'] for k, v in x['by_order'].items() if k >= 3)
lines.append('passes on which n has order 3 or more on the fold: %d' % R['passes_with_order_3_or_more'])
open('folded_alone_%s.txt' % tag, 'w').write('\n'.join(lines) + '\n'); json.dump(R, open('folded_alone_%s.json' % tag, 'w'), indent=1); print('\n'.join(lines))
