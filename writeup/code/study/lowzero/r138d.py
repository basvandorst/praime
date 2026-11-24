"""Round 138: the readings of r138_<test>.json side by side, and inverse-variance means of the amplitudes over the
fresh test Fn and the four samples of test O. usage: r138d.py"""
import os, json, math
H = os.path.dirname(os.path.abspath(__file__)); TS = ('Fn', 'R120', 'B', 'Dn', 'C'); R = {t: json.load(open(os.path.join(H, 'r138_%s.json' % t))) for t in TS}
out = {}
print('test   n    T2 p     max|m|   gamma_A          gamma_Q          gamma_off        gamma_2pt           chi2_2 chi2_3 chi2_5 (full law, 16 dof)')
for t in TS:
    r = R[t]; print('%-5s %4d  %.3f   %.5f  %.4f +- %.4f  %.4f +- %.4f  %.4f +- %.4f  %.5f +- %.5f   %5.1f %5.1f %5.1f' % (t, r['n'], r['T2_p'], max(abs(v) for v in r['m']), r['gamma_A'], r['se']['gA'], r['gamma_Q'], r['se']['gQ'],
          r['gamma_off'], r['se']['goff'], r['gamma_2pt'], r['se']['g2'], r['chi2']['2']['full'], r['chi2']['3']['full'], r['chi2']['5']['full']))
for key, sk in (('gamma_A', 'gA'), ('gamma_Q', 'gQ'), ('gamma_off', 'goff'), ('gamma_2pt', 'g2')):
    for name, ts in (('all five', TS), ('test O', TS[1:])):
        w = [1 / R[t]['se'][sk] ** 2 for t in ts]; m = sum(wi * R[t][key] for wi, t in zip(w, ts)) / sum(w); e = 1 / math.sqrt(sum(w))
        chi = sum(((R[t][key] - m) / R[t]['se'][sk]) ** 2 for t in ts); out['%s_%s' % (key, name)] = [m, e, chi]
        print('%-9s %-8s %.4f +- %.4f  (%.1f errors from 1; scatter chi2 %.1f on %d)' % (key, name, m, e, (m - 1) / e, chi, len(ts) - 1))
out['n_total'] = sum(R[t]['n'] for t in TS); print('L-functions in all: %d' % out['n_total'])
json.dump(out, open(os.path.join(H, 'r138d.json'), 'w'))
