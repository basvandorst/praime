import json, sys
import cypari2

HMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
DMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 1600
pari = cypari2.Pari()
pari.allocatemem(400 * 10**6)
out = []
for d in range(3, DMAX):
    D = -d
    if int(pari.isfundamental(D)) != 1:
        continue
    h = int(pari.qfbclassno(D))
    if h > HMAX:
        continue
    co = [int(c) for c in pari.Vec(pari.polclass(D))]
    assert co[0] == 1 and len(co) == h + 1
    out.append([D, h, [str(c) for c in co]])
out.sort(key=lambda r: (r[1], -r[0]))
json.dump(out, open('classpoly.json', 'w'))
from collections import Counter
print('discriminants', len(out), 'by class number', sorted(Counter(r[1] for r in out).items()))
print('digits in all coefficients', sum(len(c.lstrip('-')) for r in out for c in r[2]))
