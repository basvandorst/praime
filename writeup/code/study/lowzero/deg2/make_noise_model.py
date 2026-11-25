"""Noise model for the power calculation: covariance over the 20 lengths of the development residuals Sigma2_i - law_i
(20 curves, band [20, 400], mean 840 zeros), shrunk 30 % towards its diagonal; scaled per curve by 840/zeros.
usage: make_noise_model.py dev_readings.json dev_predictions.json OUT.json"""
import sys, json
import numpy as np
R = json.load(open(sys.argv[1])); P = json.load(open(sys.argv[2]))
S = np.array(R['sigma2']); F = np.array([c['full'] for c in P['curves']]); res = S - F
C = np.cov(res.T); C = 0.7 * C + 0.3 * np.diag(np.diag(C))
nz = float(np.mean(R['zeros']))
json.dump({'cov': C.tolist(), 'zeros_per_curve': nz, 'sd': np.sqrt(np.diag(C)).tolist(), 'source': 'development residuals, 20 curves, shrink 0.3'}, open(sys.argv[3], 'w'))
print('noise sd per curve, l = 0.5..10:', ' '.join('%.4f' % v for v in np.sqrt(np.diag(C))), ' mean zeros %.0f' % nz)
