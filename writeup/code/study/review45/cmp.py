import json, numpy as np
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
d = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/triple/predictions_kappa3.json'))
F = d['Fn']; idx = {D: i for i, D in enumerate(F['D'])}
my = json.load(open(SCR + 'mylaw_out.json'))
ls = np.arange(0.5, 8.01, 0.5); sel = [1, 4, 7, 10, 15]
gue = np.array(d['gue_kappa3'])
for key, r in my.items():
    D = r['D']; i = idx[D]
    full = np.array(F['full'][i]); mine = np.array(r['full'])
    print(f"D={D:6d} key={key} meanL mine {r['meanL']:.4f} json {F['meanL'][i]:.4f}")
    print('   full  json ', np.round(full[sel], 5), ' mine', np.round(mine[sel], 5), ' maxdiff(16 l) %.1e' % np.max(np.abs(full - mine)))
    print('   sigma2 maxdiff %.1e' % np.max(np.abs(np.array(F['sigma2'][i]) - r['sigma2'])),
          ' kf3 maxdiff %.1e' % np.max(np.abs(np.array(F['kf3'][i]) - r['kf3'])),
          ' kQ maxdiff %.1e' % np.max(np.abs(np.array(F['kQ'][i]) - r['kQ'])))
