"""(c) cell count vs round-13 dof; (d) linear-response estimate of the g >= 6 heuristic; (e) chi^2 of the windows."""
import json, math
import numpy as np
from scipy.stats import chi2

R = json.load(open('/home/bas/praime/writeup/data/study/r39.json'))
meas = {'train': 0.691143, 'A': 0.857295, 'B': 0.995452, 'C': 0.894136}
wins = ['train', 'A', 'B', 'C']
rdens = {w: R[w]['m'] / R[w]['n'] for w in wins}


def implied(w, mode, ref=None):
    num = den = 0.0; ncell = 0
    kap = (1 - R[w]['2']['P'][0] - R[w]['2']['per_row'][0]) / R[w]['2']['P'][0]
    for g in range(2, 39, 2):
        d = R[w][str(g)]
        for i, (P, v) in enumerate(zip(d['P'], d['per_row'])):
            ncell += 1
            if mode == 'heur' or g <= 4:
                num += v
            else:
                Pref = R[ref][str(g)]['P'][i]
                s1 = math.log(Pref / P) / math.log(rdens[ref] / rdens[w])
                num += (1 - P) - kap * P * s1 ** 2
            den += 1 - P
    P40 = [p for p in R[w]['P_pooled40'] if p > 0]
    for i, P in enumerate(R[w]['P_pooled40']):
        if P <= 0:
            continue
        ncell += 1; den += 1 - P
        if mode == 'resp40':
            Pref = R[ref]['P_pooled40'][i]
            s1 = math.log(Pref / P) / math.log(rdens[ref] / rdens[w])
            num += (1 - P) - kap * P * s1 ** 2
        else:
            num += 1 - P
    return num / den, ncell, den, kap


ref = {'train': 'A', 'A': 'train', 'B': 'A', 'C': 'B'}
print('lambda_eff', {w: round((210 / 48) / rdens[w], 2) for w in wins})
for w in wins:
    h, nc, den, kap = implied(w, 'heur')
    r1 = implied(w, 'resp', ref[w])[0]; r2 = implied(w, 'resp40', ref[w])[0]
    print(w, 'cells %d sum(1-P) %.1f kappa %.3f | implied heur %.3f (r39 %.3f) | response g>=6 %.3f | +class>=40 %.3f | measured %.3f'
          % (nc, den, kap, h, R[w]['implied_chi2_dof'], r1, r2, meas[w]))
# s-1 for a few gaps at train (ref A)
for g in (6, 8, 10, 12, 20, 30, 38):
    s = [math.log(R['A'][str(g)]['P'][i] / P) / math.log(rdens['A'] / rdens['train']) for i, P in enumerate(R['train'][str(g)]['P'])]
    print('g', g, 's-1 mean %.2f range %.2f..%.2f' % (np.mean(s), min(s), max(s)))
s = [math.log(R['A']['P_pooled40'][i] / P) / math.log(rdens['A'] / rdens['train']) for i, P in enumerate(R['train']['P_pooled40']) if P > 0]
print('class>=40 s-1 mean %.2f' % np.mean(s))
# (e) chi^2 across the four windows
pred = {w: R[w]['implied_chi2_dof'] for w in wins}
for lab, sd in (('0.058 all', {w: 0.058 for w in wins}), ('own null where known (C 0.044)', {'train': 0.058, 'A': 0.058, 'B': 0.058, 'C': 0.0443})):
    c_ind = sum(((meas[w] - 1) / sd[w]) ** 2 for w in wins)
    c_rig = sum(((meas[w] - pred[w]) / sd[w]) ** 2 for w in wins)
    print(lab, 'vs independent %.1f p=%.1e; vs rigidity %.2f p=%.3f' % (c_ind, chi2.sf(c_ind, 4), c_rig, chi2.sf(c_rig, 4)))
# noise under rigidity: each cell ~ c * chi2_1 like -> sd scales by sqrt(sum c^2 / sum (1-P)^2)
for w in wins:
    c = []; nl = []
    for g in range(2, 39, 2):
        c += R[w][str(g)]['per_row']; nl += [1 - p for p in R[w][str(g)]['P']]
    P40 = [1 - p for p in R[w]['P_pooled40'] if p > 0]; c += P40; nl += P40
    print(w, 'sd scale under rigidity %.3f' % math.sqrt(np.sum(np.square(c)) / np.sum(np.square(nl))))
sc = {'train': 0.86, 'A': 0.93, 'B': 0.95, 'C': 0.97}
c_rig = sum(((meas[w] - pred[w]) / (0.058 * sc[w])) ** 2 for w in wins)
print('vs rigidity with scaled sd: %.1f p=%.3f; train z %.2f' % (c_rig, chi2.sf(c_rig, 4), (meas['train'] - pred['train']) / (0.058 * sc['train'])))
