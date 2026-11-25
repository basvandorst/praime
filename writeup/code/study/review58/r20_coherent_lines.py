"""Predicted family-mean line amplitudes for composite r (multiplicativity: A_{r1 r2} = -A_{r1} A_{r2}/(eps c) for coprime parts)
against the measured family periodogram M(u) at u = ln r (r13_g2.npz).  Mean A_r * 150/2pi is the predicted M."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import load
B4 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4/constants4/'
B1 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind/constants/'
rows = [r for r in load() if r['cls'] == 'g2']
P = np.load('r13_g2.npz'); us, M, sd = P['us'], P['M'], P['sd']
def amp(r):
    d = json.load(open((B4 if r['rnd'] == 152 else B1) + r['key'] + '.json'))['lines']; A = {k: complex(*v) for k, v in d.items()}
    ec = r['eps'] * r['c']
    return {'2': A['2'], '1/2': A['1/2'], '3': A['3'], '1/3': A['1/3'], '6': -A['2'] * A['3'] / ec, '1/6': -A['1/2'] * A['1/3'] / ec,
            '3/2': -A['3'] * A['1/2'] / ec, '2/3': -A['2'] * A['1/3'] / ec, '1': -ec}
vals = [amp(r) for r in rows]
U = 150 / (2 * math.pi)
for k, rv in (('1', 1.0), ('2', 2.0), ('1/2', 0.5), ('3', 3.0), ('1/3', 1 / 3), ('6', 6.0), ('1/6', 1 / 6), ('3/2', 1.5), ('2/3', 2 / 3)):
    a = np.array([v[k] for v in vals]).real
    i = int(np.argmin(np.abs(us - math.log(rv))))
    # measured periodogram includes the subtraction of the main-term leakage only; at u = 0 M is the residual
    print('r = %-4s predicted family mean A_r*U = %+6.2f (sd over curves of A_r*U %.1f)   measured M(ln r) = %+6.2f%+6.2fi  (se %.2f)' % (
        k, a.mean() * U, a.std() * U, M[i].real + (a.mean() * U if k == '1' else 0), M[i].imag, sd[i]))
