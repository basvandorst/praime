"""Noise model from the USED set Fn (sigma 0.5 and 1.0): residuals of Re X and Im X after the law, decomposed into a
per-function multiplicative part (the realised one-point amplitude of each function scales its pairs), a per-function
additive part, and a pair part. The model is checked by reproducing the jackknife se of the slope on Fn; it is then
used (power.py) to estimate the se of every registered reading on the reserved designs."""
import sys, math, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for sig in (0.5, 1.0):
    S = np.load(os.path.join(HERE, 'fn_val_%.1f_P30000.npz' % sig))
    x, p = S['x'], S['predH']; n = len(S['Ds']); I, J = np.triu_indices(n, 1)
    Z = np.column_stack([np.ones(len(p)), p]); c = np.linalg.lstsq(Z, x.real, rcond=None)[0]
    rr = x.real - Z @ c; ri = x.imag - np.mean(x.imag)
    res = {}
    for nm, r in (('re', rr), ('im', ri)):
        # per-function additive effects by least squares on the pair design (r_ij = u_i + u_j + e)
        A = np.zeros((len(r), n)); A[np.arange(len(r)), I] = 1; A[np.arange(len(r)), J] = 1
        u = np.linalg.lstsq(A, r, rcond=None)[0]; e = r - A @ u
        # multiplicative: r_ij = (v_i + v_j) p_ij + ...
        B = A * p[:, None]; v = np.linalg.lstsq(np.hstack([A, B]), r, rcond=None)[0]; e2 = r - np.hstack([A, B]) @ v
        res[nm] = dict(sd=float(r.std()), sd_u=float(u.std()), sd_e=float(e.std()), sd_v=float(v[n:].std()), sd_u2=float(v[:n].std()), sd_e2=float(e2.std()))
        print('sigma %.1f %s: residual sd %.4f; additive model: sd(u) %.4f, sd(e) %.4f; add+mult: sd(u) %.4f, sd(v) %.4f, sd(e) %.4f' % (
            sig, nm, r.std(), u.std(), e.std(), v[:n].std(), v[n:].std(), e2.std()))
    # check: jackknife se of the slope under each model, simulated on the Fn design itself
    rng = np.random.default_rng(1)
    def jk(y):
        reps = []
        for f in range(n):
            k = ~((I == f) | (J == f)); reps.append(np.linalg.lstsq(Z[k], y[k], rcond=None)[0][1])
        reps = np.array(reps); return math.sqrt((n - 1) / n * ((reps - reps.mean()) ** 2).sum())
    m = res['re']
    for model in ('iid', 'add', 'addmult'):
        sl = []
        for rep in range(40):
            if model == 'iid': noise = rng.normal(0, m['sd'], len(p))
            elif model == 'add': u = rng.normal(0, m['sd_u'], n); noise = u[I] + u[J] + rng.normal(0, m['sd_e'], len(p))
            else:
                u = rng.normal(0, m['sd_u2'], n); v = rng.normal(0, m['sd_v'], n)
                noise = u[I] + u[J] + (v[I] + v[J]) * p + rng.normal(0, m['sd_e2'], len(p))
            y = p + noise; sl.append(np.linalg.lstsq(Z, y, rcond=None)[0][1])
        print('   model %-7s: sd of the slope over 40 simulations %.4f (data jackknife se: %.4f)' % (model, np.std(sl), jk(x.real)))
    out['%.1f' % sig] = res
json.dump(out, open(os.path.join(HERE, 'noise_model.json'), 'w'), indent=1)
