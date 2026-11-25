"""Synthetic zero sets for a dry run of rcross.py NN: for every discriminant of Fn, independent Poisson points with the
exact density theta'(t)/pi on [0, 600] (no orbits, no Gram lattice, no cross-correlation): every reading should be ~0.
-> synth_Fn.npz"""
import os, sys, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import theta
JS = json.load(open(os.path.join(HERE, 'predictions_cross.json'))); NEW = JS['design']['new_Fn']
rng = np.random.default_rng(7); tg = np.linspace(0.05, 600, 200001); out = {}
for D in NEW:
    q = abs(D); a = 0 if D > 0 else 1; Nt = theta(tg, q, a) / math.pi; Nt = Nt - Nt[0]
    n = rng.poisson(Nt[-1]); u = np.sort(rng.uniform(0, Nt[-1], n)); out['D%d' % D] = np.interp(u, Nt, tg)
np.savez(os.path.join(HERE, 'synth_Fn.npz'), **out); print('%d synthetic sets' % len(out))
