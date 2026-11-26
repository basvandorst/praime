# assemble predictions_quartic.json: limit predictions (pred_core.json), decisions (fixed on discovery), expected holdout errors
import sys, os, json, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); core = json.load(open(os.path.join(H, 'pred_core.json')))
ue = np.array(core['u_edges']); fi = np.round(ue * 96).astype(int); hi = slice(*core['hi_bins'])
N_HOLD_EXPECTED = {'0': 225, '1': 2150, '2': 1500}   # extrapolated from the discovery densities by halves (no holdout listing)
for key in ('0', '1', '2'):
    z = np.load(os.path.join(H, 'data', 'disc_r%s.npz' % key)); num = np.add.reduceat(z['num'], fi[:-1], axis=-1)[:, :len(ue) - 1][:, hi]
    den = np.add.reduceat(z['den'], fi[:-1], axis=-1)[:, :len(ue) - 1][:, hi]; n = len(num); N = num.sum(0); Dn = den.sum(0)
    mj = (N - num) / (Dn - den); sig = np.sqrt((n - 1) / n * ((mj - mj.mean(0))**2).sum(0))
    core['sig'][key]['sigma_disc'] = sig.tolist(); core['sig'][key]['n_disc'] = n; core['sig'][key]['n_hold_expected'] = N_HOLD_EXPECTED[key]
    core['sig'][key]['sigma_hold_expected'] = (sig * math.sqrt(n / N_HOLD_EXPECTED[key]) * 1.04).tolist()   # 1.04: ln p larger at the holdout
core['decisions'] = {'primary_sigs': ['1', '2'], 'delta_min': 25.0, 'corr_min': 0.85, 'shape_chi2_max': 43.82, 'min_fields_full_cov': 200,
                     'cond_l': [2, 3, 5, 7, 11, 13], 'cond_T': 1000.0}
core['frozen'] = '2026-11-26, before any field with 20000 <= |D| < 60000 was listed'
json.dump(core, open(os.path.join(H, 'predictions_quartic.json'), 'w'), indent=0)
import hashlib; print(hashlib.sha256(open(os.path.join(H, 'predictions_quartic.json'), 'rb').read()).hexdigest())
