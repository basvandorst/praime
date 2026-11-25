"""Verify the frozen predictions for selected triples with my own recipe code (myrecipe.py).
usage: verify20.py i0 i1   (indices into sel20.json); appends to verify20_out.jsonl"""
import sys, json, math, time, os
import numpy as np
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import myrecipe as M

P = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/moment3/predictions_triple.json'))
TT = {tuple(r['key']): r for r in P['triples']}; PP = {tuple(r['key']): r for r in P['pairs']}; SS = {r['key'][0]: r for r in P['singles']}
sel = [tuple(x) for x in json.load(open(os.path.join(HERE, 'sel20.json')))]
i0, i1 = int(sys.argv[1]), int(sys.argv[2])
out = open(os.path.join(HERE, 'verify20_out.jsonl'), 'a')
for tri in sel[i0:i1]:
    t0 = time.time()
    hm, l = M.multilinear(tri)
    hn, _ = M.multilinear(tri, no3=True)
    pf = M.P_poly(tri, hm); pn = M.P_poly(tri, hn)
    Ifull = M.int_poly_lnt(pf); Ino3 = M.int_poly_lnt(pn)
    g = M.diag_coeffs(tri)
    Idiag = quad(lambda t: M.diag_density(tri, g, t), 100, 1000, epsabs=0, epsrel=1e-11)[0]
    # pairs and singles (own)
    pairP = {}
    for a in range(3):
        for b in range(a + 1, 3):
            pr = (tri[a], tri[b]); h2, _ = M.multilinear(pr); pairP[pr] = M.P_poly(pr, h2)
    sing = {}
    for D in tri:
        h1, _ = M.multilinear((D,)); sing[D] = M.P_poly((D,), h1)
        # classical closed form for k = 1
        q = abs(D); ph = q
        for p in M.factor(q):
            ph = ph // p * (p - 1)
        cl = lambda t: ph / q * (math.log(q * t / (2 * math.pi)) + 2 * M.EG + 2 * sum(math.log(p) / (p - 1) for p in M.factor(q)))
        assert abs(sing[D](math.log(300.0)) / cl(300.0) - 1) < 1e-8, (D, sing[D](math.log(300.0)), cl(300.0))
    kirk = lambda t: np.prod([pairP[k](math.log(t)) for k in pairP]) / np.prod([sing[D](math.log(t)) for D in tri])
    Ikirk = quad(kirk, 100, 1000, epsrel=1e-11)[0]
    Ipairs = {k: M.int_poly_lnt(v) for k, v in pairP.items()}
    fr = TT[tri]
    rec = {'key': list(tri), 'full': Ifull, 'no3': Ino3, 'diag': Idiag, 'kirk': Ikirk,
           'rel_full': Ifull / fr['full'] - 1, 'rel_no3': Ino3 / fr['no3'] - 1, 'rel_diag': Idiag / fr['diag'] - 1, 'rel_kirk': Ikirk / fr['kirk'] - 1,
           'rel_pairs': [Ipairs[k] / PP[k]['full'] - 1 for k in Ipairs],
           'rel_h0': hm[()] / fr['hm'][0] - 1, 'A3ratio': hm[()] / hn[()], 'A3ratio_frozen': fr['A3_over_Aij'],
           'hm': {''.join(map(str, m)): v for m, v in hm.items()}, 'g': list(g), 'sec': time.time() - t0}
    out.write(json.dumps(rec) + '\n'); out.flush()
    print('%-22s full %+.1e no3 %+.1e diag %+.1e kirk %+.1e pairs %s h0 %+.1e A3r %.4f/%.4f (%.0f s)' % (
        tri, rec['rel_full'], rec['rel_no3'], rec['rel_diag'], rec['rel_kirk'], ' '.join('%+.1e' % x for x in rec['rel_pairs']),
        rec['rel_h0'], rec['A3ratio'], fr['A3_over_Aij'], rec['sec']), flush=True)
