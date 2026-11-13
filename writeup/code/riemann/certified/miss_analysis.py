"""Every decision of a certified chain: how hard it was, and what sits around the hard ones.
Usage: miss_analysis.py chain_decisions.json"""
import sys, json, math, collections
import numpy as np
import gmpy2
D = json.load(open(sys.argv[1]))
isp = lambda v: v > 1 and bool(gmpy2.is_prime(v))
rows = []
for d in D:
    if d['value'] is None: continue
    u = d['u']; att = d['attempts']
    ok = [a for a in att if a[1] == 'refuted']
    nodes_ok = ok[0][2] if ok else -1
    first_wrong = att[0][1] != 'refuted'
    capped = any(a[1] == 'limit' for a in att)
    rows.append(dict(u=u, j=d['j'], a=d['a'], nodes=nodes_ok, first_wrong=first_wrong, capped=capped,
                     prime=isp(u), twin=isp(u - 2) or isp(u + 2), pair_mixed=(isp(u) != isp(u + 2)) or (isp(u) != isp(u - 2))))
print(f'{len(rows)} decisions')
by_stage = collections.defaultdict(list)
for r in rows: by_stage[r['a']].append(r)
for a, rs in sorted(by_stage.items(), key=lambda kv: kv[0]):
    nd = np.array([r['nodes'] for r in rs])
    print(f'stage a = {a}: {len(rs)} decisions; nodes per proof: median {np.median(nd):.0f}, 99th percentile {np.percentile(nd, 99):.0f}, max {nd.max()}; '
          f'first LP guess wrong: {sum(r["first_wrong"] for r in rs)}; search cap hit: {sum(r["capped"] for r in rs)}')
hard_stage = max(by_stage)
rs = by_stage[hard_stage]
nd = np.array([r['nodes'] for r in rs]); thr = np.percentile(nd, 95)
hard = [r for r in rs if r['nodes'] >= thr]; easy = [r for r in rs if r['nodes'] < thr]
def frac(lst, key): return sum(r[key] for r in lst) / max(1, len(lst))
print(f'\nlast stage (a = {hard_stage}), hardest 5 % (nodes >= {thr:.0f}, {len(hard)} decisions) against the rest:')
for key in ('prime', 'twin', 'pair_mixed', 'first_wrong'):
    print(f'  {key:12s} hard {frac(hard, key):.3f}   rest {frac(easy, key):.3f}')
# position in the window
jh = np.array([r['j'] for r in hard]); je = np.array([r['j'] for r in easy])
print(f'  window position j: hard median {np.median(jh):.0f}, rest median {np.median(je):.0f}')
# local pattern: primes among u-6 .. u+6 (odd offsets)
pat = collections.Counter(''.join('1' if isp(r['u'] + o) else '0' for o in range(-6, 7, 2)) for r in hard)
print('  most common prime patterns at u-6, u-4, u-2, u, u+2, u+4, u+6 among the hardest:', pat.most_common(6))
top = sorted(rs, key=lambda r: -r['nodes'])[:40]
print('\nthe 40 hardest decisions (u: nodes, prime?, LP guess wrong?):')
print(', '.join(f"{r['u']}:{r['nodes']}{'p' if r['prime'] else ''}{'!' if r['first_wrong'] else ''}" for r in top))
seq = sorted(r['u'] for r in top[:15])
print('\nthe 15 hardest, sorted (for an OEIS lookup):', ','.join(map(str, seq)))
json.dump(dict(top=[r['u'] for r in top], seq=seq), open(sys.argv[1].replace('.json', '_analysis.json'), 'w'))
