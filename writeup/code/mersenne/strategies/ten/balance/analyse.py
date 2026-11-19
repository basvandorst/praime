import json
A = json.load(open('ref_5000.json')); B = json.load(open('ref_1e6.json')); M = json.load(open('ref_mexp_all.json'))
NAME = {'alt': 'in turn', 'light': 'on the lighter pan', 'mod4': '1 vs 3 mod 4', 'mod3': '1 vs 2 mod 3', 'mod8': '1 (x3) vs 3,5,7 mod 8', 'lastdig': 'ending 1,9 vs 3,7',
        'repeat': 'same last digit (x3) vs other', 'binsum': 'binary digit sum even vs odd', 'decsum': 'digit sum even vs odd', 'twincous': 'twins vs cousins',
        'twinsexy': 'twins (x2) vs p, p+6', 'sophie': '2p+1 vs 2p-1 prime', 'gap24': 'gap 2 vs gap 4', 'gap26': 'gap 2 vs gap 6', 'updown': 'next gap larger vs smaller',
        'liou': 'Liouville (all numbers)', 'mob': 'Mobius (square-free)', 'omega': 'distinct factors even vs odd'}
def summ(r):
    k = r['k']; d = r['L'] - r['R']
    w = 'left' if d > 1e-12 * max(1, r['L'] + r['R']) else 'right' if d < -1e-12 * max(1, r['L'] + r['R']) else 'equal'
    return dict(k=k, d=d, end=w, sL=r['aL'] / max(1, k), sR=r['aR'] / max(1, k), eq=r['eq'], ch=r['ch'], fL=r['fL'], fR=r['fR'])
def major(x): return 'left' if x['sL'] > x['sR'] else 'right' if x['sR'] > x['sL'] else 'tie'
rows = []
for key in A:
    src, rule, wt = key.split('|')
    a = summ(A[key]); b = summ(B[key]) if src == 'primes' else summ(M[key])
    held = (a['end'] == b['end']) and (major(a) == major(b))
    rows.append(dict(src=src, rule=rule, wt=wt, name=NAME[rule], a=a, b=b, held=held, same_major=major(a) == major(b), same_end=a['end'] == b['end']))
json.dump(rows, open('analysis.json', 'w'))
def pct(x): return '%d%%' % round(100 * x)
for src in ('primes', 'mexp'):
    for wt in ('one', 'val'):
        print('\n== %s, weight %s ==' % (src, wt))
        for r in rows:
            if r['src'] != src or r['wt'] != wt: continue
            a, b = r['a'], r['b']
            print('%-32s | 5k: n=%5d end %-5s L%4s R%4s eq%4d ch%4d | big: n=%7d end %-5s L%4s R%4s eq%6d ch%6d | %s' % (r['name'], a['k'], a['end'], pct(a['sL']), pct(a['sR']), a['eq'], a['ch'], b['k'], b['end'], pct(b['sL']), pct(b['sR']), b['eq'], b['ch'], 'holds' if r['held'] else ('lead flips' if not r['same_end'] else 'majority flips')))
held = sum(r['held'] for r in rows); print('\nconclusions (winner at the end and majority of the time) that held:', held, 'of', len(rows))
