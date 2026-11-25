"""Swapped-term predictions for the near-equal-conductor pairs of a design (theory only).
usage: swap_run.py DESIGN START STOP   (DESIGN in NN NO OO; pairs with averaging factor F > 0.1, sorted by (D1, D2))
-> swap_DESIGN_START.npz with D1, D2, F, C (len x 16)"""
import os, sys, json, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
from swaplib import cswap
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..')); LOW = os.path.join(REPO, 'writeup/code/study/lowzero')
PJ = json.load(open(os.path.join(LOW, 'triple', 'predictions_kappa3.json')))
OLD = PJ['R120']['D'] + PJ['B']['D'] + PJ['Dn']['D']; NEW = PJ['Fn']['D']
T1, T2, H = 100.0, 600.0, 0.02; W = 0.3 * np.arange(1, 17)
PMAX = int(os.environ.get('SWAP_PMAX', '200000')); DU = float(os.environ.get('SWAP_DU', '0.03')); PH = primes_upto(PMAX)
def Fapprox(a, b):
    d = math.log(abs(b) / abs(a)); x = 250 * d
    return 1.0 if x == 0 else abs(math.sin(x) / x)
def pairs(design):
    out = []
    if design == 'NN':
        for i in range(len(NEW)):
            for j in range(i + 1, len(NEW)): out.append((NEW[i], NEW[j]))
    elif design == 'NO':
        out = [(a, b) for a in NEW for b in OLD]
    else:
        for i in range(len(OLD)):
            for j in range(i + 1, len(OLD)): out.append((OLD[i], OLD[j]))
    return [(a, b, Fapprox(a, b)) for a, b in out if abs(a) != abs(b) and Fapprox(a, b) > 0.1]
if __name__ == '__main__':
    des, s0, s1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    PL = pairs(des); print('%s: %d pairs with F > 0.1 (equal |D| excluded: their term is exactly 0)' % (des, len(PL)), flush=True)
    sub = PL[s0:s1]; X = {}
    Ds = sorted(set([a for a, b, f in sub] + [b for a, b, f in sub])); XT = kron_table(Ds, PH); X = {D: XT[i] for i, D in enumerate(Ds)}
    C = np.zeros((len(sub), len(W))); t0 = time.time()
    for n, (a, b, f) in enumerate(sub):
        C[n] = cswap(a, b, W, T1, T2, H, PH, X[a], X[b], du=DU)
    np.savez(os.path.join(HERE, 'swap_%s_P%d_%05d.npz' % (des, PMAX, s0)), D1=np.array([a for a, b, f in sub]), D2=np.array([b for a, b, f in sub]), F=np.array([f for a, b, f in sub]), C=C)
    print('pairs %d-%d done in %.0fs' % (s0, s0 + len(sub), time.time() - t0))
