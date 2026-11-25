import json, math, sys
import numpy as np
D = json.load(open('all_decomp.json'))
print('%-16s %5s | %-17s | %-13s | %-17s | %s' % ('window', 'Q', 'u2X meas/pred', 'u2Z meas/pred', 'u2R meas/pred', 'beta'))
for w, d in D.items():
    for Q in (210, 2310, 30030):
        e = d['h2_%d' % Q]
        print('%-16s %5d | %.3f/%.3f %.2f | %.3f/%.3f | %.3f/%.3f %.2f | %.2f+-%.2f' % (w, Q, e['measured_u2X'], e['pred_u2X'], e['ratio_X'],
              e['measured_u2Z'], e['pred_u2Z'], e['measured_u2R'], e['pred_u2R'], e['ratio_R'], e['beta_hat'], e['se_beta_approx']))
