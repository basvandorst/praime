# compute the predicted bin averages (two truncation levels: the difference is the numerical error) and the comparators
import os, sys, time, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import predlib as PL
edges = np.round(np.linspace(0, 3, 61), 10); out = {'edges': edges.tolist()}
for sgn in (-1, 1):
    t0 = time.time()
    m1, tot1 = PL.predict_bins(edges, sgn, Amax=1e6, tau=1e-6, B0=500, Bcap=10000, verbose=True)
    m0, tot0 = PL.predict_bins(edges, sgn, Amax=3e5, tau=3e-6, B0=300, Bcap=5000)
    k = '%+d' % sgn
    out[k] = {'pred': m1.tolist(), 'numerr': np.abs(m1 - m0).tolist(), 'lines_no_square_term': tot1.tolist(),
              'single_line_zeta3': PL.single_line_bins(edges, sgn).tolist(), 'single_line_unit': PL.single_line_bins(edges, sgn, 1.0).tolist(),
              'gl2_prime_level': PL.gl2_prime_level_bins(edges, sgn).tolist(), 'tt_secondary': PL.tt_bins(edges, sgn).tolist()}
    print('sgn %s: %.1f s; max numerical difference %.4f, rms %.4f' % (k, time.time() - t0, np.max(np.abs(m1 - m0)), np.sqrt(np.mean((m1 - m0)**2))))
json.dump(out, open('pred_cache.json', 'w'))
