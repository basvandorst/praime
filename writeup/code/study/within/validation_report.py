"""Prints the validation-stretch measurements (data/val_measure.json, data/val_tail.json; stretch [1.5e11, 1.5e11+1e8)
only) next to the model components (data/comp_25.7342.json and the y = 13 tail bins). No within part is assembled."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json, math
from tailmod import responses_all, autocorr, tail_values
vm = json.load(open('data/val_measure.json')); vt = json.load(open('data/val_tail.json'))
lam = vm['lam']; print('validation stretch, lam = %.4f, %d primes' % (lam, vm['nprimes']))
b = vm['b']
print('q = 3: b_j = E[chi(p_n) chi(p_{n+j})]:  ' + '  '.join('%s:%+.5f(%.5f)' % (j, b[j][0], b[j][1]) for j in list(b)[:10]))
print('       LOS: b_1 = -ln(2 pi lam/3)/(2 lam) = %+.5f, b_j = -1/(2(j-1) lam): %s' % (-math.log(2 * math.pi * lam / 3) / (2 * lam), ' '.join('%+.5f' % (-1 / (2 * (j - 1) * lam)) for j in range(2, 8))))
print('       stays rho_j (k=1): ' + ' '.join('%s:%+.5f' % (j, v[0]) for j, v in list(vm['rho_k1'].items())[:6]) + '   identity (b2-b1^2)/(1-b1^2) = %+.5f' % ((b['2'][0] - b['1'][0] ** 2) / (1 - b['1'][0] ** 2)))
print('       connected 4-point parts (k=1, beyond pairings): ' + ' '.join('%s:%+.5f(%.5f)' % (j, v['conn'], v['E4'][1]) for j, v in list(vm['c4_k1'].items())[:5]))
BINS = [tuple(x) for x in vt['bins']]
Q, U, R = responses_all(13, lam, (3, 5, 7))
for q in (3, 5, 7):
    A = autocorr(Q, U, R[q]['r']); tv = tail_values(Q, U, A, R[q]['C'], lam, 13, [10 ** 4], BINS)
    for k in (1, 2):
        mb = vt['%d_%d' % (q, k)]['binsum']
        line = '  q=%d k=%d lag-bin sums of tau_j/(C-1)^2:' % (q, k)
        line += '  '.join(' [%d,%d) %+.5f(%.5f)' % (lo, hi, m[0], m[1]) + (' tail %+.5f' % tv['bins'][(lo, hi)] if k == 1 and lo >= 11 else '') for (lo, hi), m in zip(BINS, mb))
        print(line)
