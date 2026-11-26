# frozen predictions for the S4 quartic murmuration test: limit density on the registered bins, comparators, numerical error
import sys, os, json, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P
H = os.path.dirname(os.path.abspath(__file__))
ue = np.array([0.0] + [(2 * j + 1) / 24 for j in range(24)]); ye = ue**3
t0 = time.time(); M = {T: P.quartic_model(T=T) for T in (750.0, 1500.0)}; cub = P.cubic_model(T=1500.0); one = P.single_line_model(P.A1_exact(), T=1500.0)
out = {'u_edges': ue.tolist(), 'y_edges': ye.tolist(), 'hi_bins': [4, 24], 'A1': P.A1_exact(),
       'bhargava_c': {'0': 1 / 48 / P.A1_exact(), '1': 1 / 8 / P.A1_exact(), '2': 1 / 16 / P.A1_exact()}, 'sig': {}}
for r2, mus in P.SIG_MU.items():
    sg = -1 if r2 == 1 else 1; p15 = P.bin_means(M[1500.0], mus, ye); p7 = P.bin_means(M[750.0], mus, ye)
    out['sig'][str(r2)] = {'mu': list(mus), 'pred': p15.tolist(), 'numerr': np.abs(p15 - p7).tolist(), 'tail_y^-1/3': P.tail_coeffs(M[1500.0], mus)['0.3333'],
        'comparators': {'one_line': P.bin_means(one, mus, ye, squares=False).tolist(), 'deg2_kernel': P.bin_means(M[1500.0], P.CUBIC_MU[sg], ye).tolist(),
                        'cubic_transplant': P.bin_means(cub, P.CUBIC_MU[sg], ye).tolist(), 'no_square_term': P.bin_means(M[1500.0], mus, ye, squares=False).tolist()}}
    print('r2=%d: max numerical error (T 1500 vs 750) on bins 1..23: %.1e' % (r2, np.max(np.abs(p15 - p7)[1:])), flush=True)
out['build_seconds'] = time.time() - t0
json.dump(out, open(os.path.join(H, 'pred_core.json'), 'w'), indent=0)
