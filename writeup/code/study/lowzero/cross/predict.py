"""Frozen predictions for the cross-correlation test (round_draft.txt). Uses only the discriminants (no zeros).
Law: diagonal (non-swapped) ratios term = semiclassical shared-orbit sum, in its finite-range form on the registered grid
(orbits n <= 300 with the exact grid averages of the cosine products; orbits n > 300 in the T -> infinity form, primes
<= 2e5). Writes predictions_cross.json (tables, expected readings under the hypotheses) and predictions_cross_pairs.npz
(per-pair FULL, REST, psi at 2,3,5,7 for the designs NN = Fn x Fn and NO = Fn x old), prints sha256 of both.
Statistic: w_k = 0.3 k (k = 1..16); t-grid T1 + 0.02 j, t <= 600 - w; dn_j(t) = #{gamma_j in (t, t+w]} - (theta_j(t+w)
- theta_j(t))/pi; C_ij(w) = mean_t dn_i dn_j (no mean subtraction)."""
import os, sys, json, math, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
from readings import readings
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..')); LOW = os.path.join(REPO, 'writeup/code/study/lowzero')
PJ = json.load(open(os.path.join(LOW, 'triple', 'predictions_kappa3.json')))
OLD = PJ['R120']['D'] + PJ['B']['D'] + PJ['Dn']['D']; NEW = PJ['Fn']['D']; ALL = NEW + OLD; nN = len(NEW)
T1, T2, H = 100.0, 600.0, 0.02; W = 0.3 * np.arange(1, 17); PR = (2, 3, 5, 7); X0 = 300
def grid(w):
    M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1; return T1 + H * np.arange(M)
def AT(p, w):
    tg = grid(w); lp = math.log(p); A = 0.0; B = 0.0; A0 = 0.0; B0 = 0.0; m = 1
    while 1.0 / (m * m * p ** m) > 1e-14:
        base = 2 / math.pi ** 2 / (m * m * p ** m) * math.sin(m * w * lp / 2) ** 2
        f = 1 + float(np.mean(np.cos(2 * m * (tg + w / 2) * lp)))
        if m % 2: A += base * f; A0 += base
        else: B += base * f; B0 += base
        m += 1
    return A, B, A0, B0
TAB = np.array([[AT(p, w) for p in PR] for w in W])          # (16, 4, 4): A_T, B_T, A_inf, B_inf
A_T = TAB[:, :, 0]; B_T = TAB[:, :, 1]
P = primes_upto(200000); XA = kron_table(ALL, P)
Cinf = pair_cov_full(XA, W, P); REST = pair_cov_full(XA[:, 4:], W, P[4:])
# finite-range correction for orbits n <= X0
PP = [int(p) for p in P if p <= X0]; PW = []
for j, p in enumerate(PP):
    m = 1
    while p ** m <= X0: PW.append((p ** m, m, j)); m += 1
lnn = np.array([math.log(n) for n, _, _ in PW]); mm = np.array([m for _, m, _ in PW]); jj = np.array([j for _, _, j in PW])
CHIN = XA[:, :len(PP)].astype(float)[:, jj] ** mm[None, :]
DELTA = np.zeros_like(Cinf)
for k, w in enumerate(W):
    tg = grid(w); Cs = np.cos(np.outer(tg + w / 2, lnn)); K = Cs.T @ Cs / len(tg)
    c = -(2 / math.pi) / (mm * np.exp(lnn / 2)) * np.sin(w * lnn / 2); V = CHIN * c[None, :]
    DELTA[k] = V @ K @ V.T - (CHIN * (c ** 2 / 2)[None, :]) @ CHIN.T
FULL = Cinf + DELTA
designs = {'NN': np.triu_indices(nN, 1)}
I2, J2 = np.meshgrid(np.arange(nN), np.arange(nN, len(ALL)), indexing='ij'); designs['NO'] = (I2.ravel(), J2.ravel())
I3, J3 = np.triu_indices(len(OLD), 1); designs['OO'] = (I3 + nN, J3 + nN)
out = {}; expected = {}; rng = np.random.default_rng(140)
for name, (I, J) in designs.items():
    F = np.array([FULL[k][I, J] for k in range(len(W))]); R = np.array([REST[k][I, J] for k in range(len(W))])
    PSI = np.column_stack([XA[I, a].astype(float) * XA[J, a].astype(float) for a in range(4)])
    habs = np.array([sum((A_T[k, a] + B_T[k, a]) * PSI[:, a] ** 2 for a in range(4)) for k in range(len(W))])
    e = {'H_law': readings(F, PSI, R, F, A_T, B_T), 'H_indep': readings(0 * F, PSI, R, F, A_T, B_T),
         'H_flip': readings(-F, PSI, R, F, A_T, B_T), 'H_abs': readings(habs, PSI, R, F, A_T, B_T),
         'H_law_Tinf': readings(np.array([Cinf[k][I, J] for k in range(len(W))]), PSI, R, F, A_T, B_T),
         'dry_run_law_plus_noise_0.003': readings(F + 0.003 * rng.standard_normal(F.shape), PSI, R, F, A_T, B_T)}
    for h in e:
        keep = (h == 'H_law')
        for key in ('icpt', 's', 'u'):
            v = e[h].pop(key, None)
            if keep: e[h]['per_w_' + key] = v
    expected[name] = e
    print('%s: %d pairs; sd over pairs of the law per w: %s' % (name, len(I), ' '.join('%.4f' % s for s in F.std(1))))
    print('  finite-range part: sd over pairs per w: %s' % ' '.join('%.4f' % s for s in np.array([DELTA[k][I, J].std() for k in range(len(W))])))
    for h, r in e.items(): print('  expected readings %-28s ' % h + ', '.join('%s %+.4f' % (k, v) for k, v in r.items() if not k.startswith('per_w')))
    if name != 'OO':
        out[name + '_full'] = F.T.astype(np.float32); out[name + '_rest'] = R.T.astype(np.float32); out[name + '_psi'] = PSI.astype(np.int8)
        out[name + '_D1'] = np.array(ALL)[I]; out[name + '_D2'] = np.array(ALL)[J]
fn = os.path.join(HERE, 'predictions_cross_pairs.npz'); np.savez_compressed(fn, **out)
sha_npz = hashlib.sha256(open(fn, 'rb').read()).hexdigest()
JS = {'statistic': 'C_ij(w) = mean over t in T1 + 0.02 j (t <= 600 - w) of dn_i(t) dn_j(t); dn = #zeros in (t, t+w] - (theta(t+w) - theta(t))/pi, theta(t) = (t/2) ln(q/pi) + Im lnGamma((1/2 + a + it)/2)',
      'T1': T1, 'T2': T2, 'h': H, 'w': W.tolist(), 'primes': PR, 'X0_finite_range': X0, 'pmax_pairs': 200000,
      'A_T': A_T.T.tolist(), 'B_T': B_T.T.tolist(), 'A_inf': TAB[:, :, 2].T.tolist(), 'B_inf': TAB[:, :, 3].T.tolist(),
      'design': {'new_Fn': NEW, 'old': OLD, 'old_sets': 'R120 (r117/r120 parts) + B (r134_B) + Dn (r136_Dn)'},
      'pairs_file': 'predictions_cross_pairs.npz', 'pairs_sha256': sha_npz, 'expected_readings': expected,
      'hypotheses': {'H_law': 'shared-orbit diagonal law, finite-range form: all amplitudes 1', 'H_indep': 'no cross-correlation: all 0',
                     'H_flip': 'the diagonal term with the wrong sign: all -1', 'H_abs': 'shared orbits without their phases: Cov = sum_p (A_p + B_p) 1[p not | k1k2]',
                     'H_law_Tinf': 'the law without the finite-range terms (to show their size)'}}
json.dump(JS, open(os.path.join(HERE, 'predictions_cross.json'), 'w'))
sha_json = hashlib.sha256(open(os.path.join(HERE, 'predictions_cross.json'), 'rb').read()).hexdigest()
print('predictions_cross_pairs.npz sha256 %s' % sha_npz); print('predictions_cross.json sha256 %s' % sha_json)
for a, p in enumerate(PR):
    print('p=%d A_T ' % p + ' '.join('%.5f' % v for v in A_T[:, a])); print('    B_T ' + ' '.join('%.5f' % v for v in B_T[:, a]))
    print('    A_T/A_inf - 1 ' + ' '.join('%+.4f' % (x / y - 1) for x, y in zip(A_T[:, a], TAB[:, a, 2])))
