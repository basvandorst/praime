"""Order-3 + overlap correction delta(lam) = ln(share3/share2) at several lam, against r53.json's three values, and the
effect of interpolating delta linearly on [30, 38, 42.84] versus a smooth fit through more points."""
import sys, math, json, time
import numpy as np
from law import Gap

r53 = json.load(open('/home/bas/praime/writeup/data/study/r53.json'))
data = json.load(open('/home/bas/praime/writeup/data/study/oes_gap_counts_4e18.json'))['N']
L3 = np.array([30.0, 38.0, 42.84])
B = math.log(4e18)
lams = [float(x) for x in sys.argv[2].split(',')]
out = {}
for g in map(int, sys.argv[1].split(',')):
    t0 = time.time()
    G = Gap(g)
    d = G.delta3(lams)
    # weighted mean of delta over the integrand (weight share2 * e^u / u) with lam in [30, B], fine u grid:
    u = np.linspace(30, B, 801)
    ls2 = np.array([G.lnshare2(x) for x in u])
    wgt = np.exp(ls2 + u - ls2.max() - B) / u
    lin_r53 = np.interp(u, L3, r53['o3'][str(g)]) if str(g) in r53['o3'] else None
    lin_own3 = np.interp(u, L3, [d[lams.index(x)] for x in L3]) if all(x in lams for x in L3) else None
    # smooth: polynomial in 1/lam through all computed points (degree len-1, at most 3)
    deg = min(3, len(lams) - 1)
    cf = np.polyfit(1 / np.array(lams), d, deg)
    smooth = np.polyval(cf, 1 / u)
    W = wgt.sum()
    row = {'lams': lams, 'delta_own': d, 'delta_r53': r53['o3'].get(str(g)),
           'wmean_smooth': float((wgt * smooth).sum() / W),
           'wmean_lin_own3': None if lin_own3 is None else float((wgt * lin_own3).sum() / W),
           'wmean_lin_r53': None if lin_r53 is None else float((wgt * lin_r53).sum() / W),
           'secs': time.time() - t0}
    out[g] = row
    print(g, json.dumps(row), flush=True)
json.dump(out, open('/home/bas/scratch/study-review12/o3_%s.json' % sys.argv[1].replace(',', '_'), 'w'), indent=1)
