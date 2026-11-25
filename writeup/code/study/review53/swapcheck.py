"""Direct 8-term swap sum at shifts z = eps (1, 1.31, 1.73), eps -> 0 by polynomial extrapolation, vs the multilinear
formula (2.2) with my coefficients (verify20_out.jsonl) and the frozen P_at values."""
import sys, os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import myrecipe as M
P = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/moment3/predictions_triple.json'))
TT = {tuple(r['key']): r for r in P['triples']}
mine = {tuple(json.loads(l)['key']): json.loads(l) for l in open(os.path.join(HERE, 'verify20_out.jsonl'))}
for tri in [tuple(int(a) for a in s.split(',')) for s in sys.argv[1:]]:
    hm = {tuple(int(c) for c in k): v for k, v in mine[tri]['hm'].items()}
    epss = [0.004, 0.002, 0.001, 0.0005]
    for t in (100.0, 1000.0):
        vals = [M.swap_sum(tri, t, e) for e in epss]
        co = np.polyfit(epss, vals, len(epss) - 1)
        ext = co[-1]
        mlP = M.P_at(tri, hm, t)
        fr = TT[tri]['P_at']['full'][TT[tri]['P_at']['t'].index(t)]
        print('%s t %5.0f swap sums %s -> %.6f | (2.2) mine %.6f frozen %.6f | rel %.1e %.1e' % (
            tri, t, np.round(vals, 3).tolist(), ext, mlP, fr, ext / mlP - 1, fr / mlP - 1), flush=True)
