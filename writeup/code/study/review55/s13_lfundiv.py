import os, json, math, time, sys
import numpy as np, cypari2
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'
P = json.load(open(B + '/predictions_blind.json')); zz = np.load(sys.argv[1])
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
for key in sys.argv[2:]:
    e = P['products'][key]; t0 = time.time()
    L = pari('lfundiv(lfuncreate(%s), lfuncreate(1))' % e['spec'][1])
    Li = pari.lfuninit(L, [1, 303])
    gd = np.sort(np.array([float(x) for x in pari.lfunzeros(Li, [0.001, 300.0], 16)]))
    g = np.load(B + '/zeros_blind/' + key + '.npz')['g']
    isz = np.array([np.min(np.abs(zz - v)) < 1e-6 for v in g]); gr = g[~isz]; gr = gr[(gr > 0) & (gr <= 300)]
    same = len(gd) == len(gr) and np.max(np.abs(gd - gr)) < 1e-6
    print('%s: direct lfundiv zeros (divz 16) %d, subtraction list %d, identical to 1e-6: %s, max diff %.1e, feq %s, %.0fs' % (
        key, len(gd), len(gr), same, np.max(np.abs(gd - gr)) if len(gd) == len(gr) else float('nan'), pari.lfuncheckfeq(L), time.time() - t0), flush=True)
