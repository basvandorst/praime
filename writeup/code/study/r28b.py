"""Round 28 diagnostic (after the round-27 data were seen): the train tails >= 6, 8, 10 ln x with the order-4 terms of
round 8 (spans 122..200) replacing order 3, plus the overlap correction, against the counts of round 27."""
import os, json, math, pickle
import numpy as np
import theory2 as T2
import theory_r3 as R3
import r8

OUT = R3.OUT
c3 = pickle.load(open(R3.CACHE, 'rb')); c8 = pickle.load(open(r8.CACHE8, 'rb'))
L, W = R3.WIN['train']
u = np.linspace(math.log(L), math.log(L + W), 401); w = np.exp(u) / u
r27 = json.load(open(os.path.join(OUT, 'r27.json')))
out = {}
for variant in ('order3_beyond_120', 'order4_to_200'):
    gl = {}
    for h in range(2, R3.HMAX + 1, 2):
        SB, la = c3[(0, h)]
        if SB == 0:
            gl[h] = 0.0; continue
        lav = np.interp(u, R3.LAMG, la)
        pt = T2.Pattern([0, h], 2)
        lav = lav + np.interp(u, R3.LAMG, [T2.overlap_correction(pt, l) for l in R3.LAMG])
        if variant == 'order4_to_200' and 120 < h <= 200 and c8.get(h) is not None:
            lav = lav + np.interp(u, r8.LAMT, c8[h][0])
        gl[h] = float(np.trapezoid(w * SB / u * np.exp(lav), u) / np.trapezoid(w, u))
    tot = sum(gl.values())
    out[variant] = {}
    for k in (6, 8, 10):
        d = r27['train_%dlnx' % k]
        law_share = sum(v for h, v in gl.items() if h >= d['threshold']) / tot
        n = d['expected'] / (sum(v for h, v in gl.items() if h >= d['threshold']) / tot) if variant == 'x' else None
        out[variant]['%dlnx' % k] = law_share
print(json.dumps(out, indent=1))
for k in (6, 8, 10):
    d = r27['train_%dlnx' % k]; base = out['order3_beyond_120']['%dlnx' % k]; new = out['order4_to_200']['%dlnx' % k]
    print(k, 'count', d['count'], 'expected before %.1f after %.1f' % (d['expected'], d['expected'] * new / base), 'ratio after %.3f' % (d['count'] / (d['expected'] * new / base)))
json.dump(out, open(os.path.join(OUT, 'r28b.json'), 'w'), indent=1)
