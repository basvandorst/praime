import json, numpy as np
for lf, cf, jf in [('law1e9.npy','c1e9.npz','r85_1e9.json'),('law1e12.npy','c1e12.npz','r85_1e12.json'),('law1e15.npy','c1e15.npz','r86_1e15.json')]:
    T = np.load(lf); r = json.load(open('/home/bas/praime/writeup/data/study/' + jf)); cnt = np.load(cf)['cnt']; NB = cnt.shape[1]
    zb_all = []; zj_all = []; zP = []; zU = []
    for k in (1, 2, 3):
        law = T[k-1] / T[k-1].sum(1, keepdims=True); J = np.array(r['k'][str(k)]['law'])
        tot = cnt[k-1].sum(0).reshape(4, 4).astype(float); N = tot.sum(1)[:, None]; dat = tot / N
        sebin = np.sqrt(dat * (1 - dat) / N)
        bl = cnt[k-1].reshape(NB, 4, 4).astype(float); bp = bl / bl.sum(2)[:, :, None]
        sebm = bp.std(0, ddof=1) / np.sqrt(NB)
        zb = (law - dat) / sebin; zm = (law - dat) / sebm
        zb_all += list(zb.ravel()); zj_all += list(np.array(r['k'][str(k)]['z']).ravel())
        # Pearson multinomial chi2 per row (3 dof per row)
        zP.append(float((N * (dat - law) ** 2 / law).sum()))
        zU.append(float(np.abs((0.25 - dat) / sebin).max()))
        print(jf, 'k', k, 'max|law ours - law study| %.1e' % np.abs(law - J).max(), 'max|z| bin %.2f batch %.2f' % (np.abs(zb).max(), np.abs(zm).max()),
              'Pearson chi2 %.1f (12 dof)' % zP[-1], 'uniform max|z| %.0f' % zU[-1])
    print('   all: max|z| %.3f (json %.3f), sum z^2 %.2f (json %.2f), Pearson total %.1f on 36 dof' % (np.abs(zb_all).max(), np.abs(zj_all).max(), np.sum(np.square(zb_all)), np.sum(np.square(zj_all)), sum(zP)))
