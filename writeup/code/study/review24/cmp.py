import sys, json, numpy as np
c = np.load(sys.argv[1]); r = json.load(open(sys.argv[2]))
cnt = c['cnt']; NB = cnt.shape[1]
for k in (1, 2, 3):
    tot = cnt[k-1].sum(0).reshape(4, 4).astype(float)
    N = tot.sum(1); dat = tot / N[:, None]
    J = r['k'][str(k)]
    print('k', k, 'rows equal:', np.allclose(N, J['row_counts'], rtol=0, atol=0), N, 'max|dat diff| %.2e' % np.abs(dat - np.array(J['data'])).max())
    # batch means SE
    bl = cnt[k-1].reshape(NB, 4, 4).astype(float)
    bn = bl.sum(2); bp = bl / bn[:, :, None]
    se_bm = bp.std(0, ddof=1) / np.sqrt(NB)
    se_bin = np.sqrt(dat * (1 - dat) / N[:, None])
    print('   ratio batch-means SE / binomial SE: mean %.3f  min %.3f max %.3f' % ((se_bm/se_bin).mean(), (se_bm/se_bin).min(), (se_bm/se_bin).max()))
    print('   ratio^2 per row:', np.round(((se_bm/se_bin)**2).mean(1), 3))
