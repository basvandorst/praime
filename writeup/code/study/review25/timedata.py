import sys, time; sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import r89
t0 = time.time(); c = r89.data(); print('data', round(time.time() - t0, 1), 's', {k: sum(v.values()) for k, v in c.items()})
import json; mine = json.load(open('counts.json'))
for k in (1, 2, 3): assert {int(a): b for a, b in mine[f'{k}_in'].items()} == c[k], k
print('per-distance counts identical to own recount')
