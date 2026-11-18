"""Timings of the test with a receipt, each the smaller of two runs (the machine is shared), with the exact counts of operations."""
import json, sys, time
from receipt import One, Two, run, build, check

rows = []
lines = []
for p in [int(x) for x in sys.argv[1].split(',')]:
    for cls, name in ((One, 'one cell'), (Two, 'two cells')):
        for k in (6, 8):
            best = None
            for rep in range(2):
                G = cls(p)
                t = time.time(); y, saved, top = run(G, k); t1 = time.time() - t
                cnt = [0]
                t = time.time(); mus = build(G, k, y, saved, count=cnt); t2 = time.time() - t
                c2 = [0]
                t = time.time(); ok = check(G, y, mus, count=c2); t3 = time.time() - t
                assert ok and G.verdict(y)
                cur = [t1, t2, t3]
                best = cur if best is None else [min(a, b) for a, b in zip(best, cur)]
            row = dict(p=p, group=name, k=k, squarings=top, saved=len(saved), receipt_products=cnt[0], check_squarings=c2[0], small_powers=2 * k, test_s=round(best[0], 3), receipt_s=round(best[1], 3), check_s=round(best[2], 4), ratio=round(best[0] / best[2], 1))
            rows.append(row)
            s = 'p = %d, %s, %d halvings: test %d squarings %.2f s; receipt from %d saved residues %.2f s; check %d squarings and %d small powers %.3f s; ratio %.0f' % (p, name, k, top, best[0], len(saved), best[1], c2[0], 2 * k, best[2], best[0] / best[2])
            print(s, flush=True)
            lines.append(s)
json.dump(rows, open('receipt_timing.json', 'w'))
open('receipt_timing.txt', 'w').write('\n'.join(lines) + '\n')
