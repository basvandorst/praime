"""Second-order check of the lag-1 term: LOS's first-order non-Markov tensor on top of the actual transition matrix
(model B):  P(c | a, b) = P_bc [1 + (phi/(2 lam)) (P_ba - delta_ac)],  D3(a,b,c) = pi_a P_ab P(c|a,b).
P taken from the validation stretch (a proxy for the HL P, which agrees with it: B_q within 0.0006)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, json, math
from triple import tau1_from_D3
vm = json.load(open('data/val_measure.json')); lam = vm['lam']
for q in (3, 5, 7):
    P = np.array(vm['tables']['%d_1' % q]['P']); pi = np.array(vm['tables']['%d_1' % q]['pi']); C = len(P)
    D3 = np.zeros((C, C, C))
    for a in range(C):
        for b in range(C):
            for c in range(C):
                D3[a, b, c] = pi[a] * P[a, b] * P[b, c] * (1 + C / (2 * lam) * (P[b, a] - (a == c)))
    print('q=%d  model B tau1/(C-1)^2 = %+.5f   leading -1/(2(C-1)lam) = %+.5f   (HL triple law -> see law_output)' % (q, tau1_from_D3(D3)[0], -1 / (2 * (C - 1) * lam)))
