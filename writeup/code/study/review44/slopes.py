import json, math, numpy as np
A = json.load(open('agg.json')); TE = json.load(open('tailemp.json'))
P = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/within/data/predictions.json'))
OLDF = [json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r135b_%s.json' % n)) for n in ('freshA', 'reg', 'freshB')]
cmpn = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/within/data/comp_22.6973.json'))
cmpo = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/within/data/comp_27.1203.json'))
HN = [10**4, 2*10**4, 5*10**4, 10**5, 2*10**5, 5*10**5, 10**6]; HO = [10**4, 10**5, 10**6]
def fit(Wmat, Hs):
    """Wmat: windows x H. weighted LS slope of the pooled means vs ln H (weights from the per-H window scatter); jackknife error"""
    x = np.log(np.array(Hs, float)); n = Wmat.shape[0]
    def slope(M):
        y = M.mean(0); s2 = M.var(0, ddof=1)/M.shape[0]; w = 1/s2
        xb = (w*x).sum()/w.sum(); return ((w*(x-xb)*(y-(w*y).sum()/w.sum())).sum())/((w*(x-xb)**2).sum())
    s = slope(Wmat); jk = np.array([slope(np.delete(Wmat, i, 0)) for i in range(n)])
    return s, math.sqrt((n-1)/n*((jk-jk.mean())**2).sum())
def diffs(Wmat, Hs, a, b):
    d = (Wmat[:, Hs.index(b)] - Wmat[:, Hs.index(a)])/math.log(b/a); return d.mean(), d.std(ddof=1)/math.sqrt(len(d))
rows = {}
for q, k in ((3,1),(5,1),(7,1),(3,2),(10,2)):
    qq = 5 if q == 10 else q; C = len([a for a in range(qq) if math.gcd(a, qq) == 1])
    Wn = np.array([A['win_%d_%d_%d' % (q, k, H)]['pw'] for H in HN]).T
    Wo = np.array([sum((J['win_%d_%d_%d' % (q, k, H)]['per_window'] for J in OLDF), []) for H in HO]).T
    sn, en = fit(Wn, HN); so, eo = fit(Wo, HO)
    lawn = (P['%d_%d_%d' % (qq, k, 10**6)]['law'] - P['%d_%d_%d' % (qq, k, 10**4)]['law'])/math.log(100)
    lawo = (P['%d_%d_%d' % (qq, k, 10**6)]['old_law'] - P['%d_%d_%d' % (qq, k, 10**4)]['old_law'])/math.log(100)
    Vn = cmpn[str(qq)]['V_k%d' % k]; Vo = cmpo[str(qq)]['V_k%d' % k]
    s = 'q=%2d k=%d  NEW: fit slope/unit lnH %+.5f +- %.5f   law 1e4-1e6 %+.5f (-V/lam %+.5f)  ratio %.2f  z %+.1f' % (q, k, sn, en, lawn, -Vn/cmpn['lam'], sn/lawn, (sn-lawn)/en)
    if k == 1:
        te = TE['%d_17' % q]['W_tail']; emp17 = (te['1000000'] - te['10000'])/math.log(100)
        s += '  | empirical-response y=17 tail slope %+.5f  ratio %.2f  z %+.1f' % (emp17, sn/emp17, (sn-emp17)/en)
    print(s)
    print('          diffs per unit lnH new: ' + '  '.join('[%g,%g] %+.5f+-%.5f' % ((a, b) + diffs(Wn, HN, a, b)) for a, b in ((10**4, 10**5), (10**5, 10**6), (10**4, 10**6))))
    print('          OLD (lam 27.0): fit slope %+.5f +- %.5f   law %+.5f  ratio %.2f  z %+.1f;  new/old measured slope ratio %.2f +- %.2f (1/lam^2 predicts %.2f)'
          % (so, eo, lawo, so/lawo, (so-lawo)/eo, sn/so, abs(sn/so)*math.hypot(en/sn, eo/so), (27.017/22.693)**2))
    rows['%d_%d' % (q, k)] = dict(new=(sn, en), old=(so, eo), law_new=lawn, law_old=lawo)
# pooled q = 5, 7 (k = 1)
for which in ('new', 'old'):
    s = np.mean([rows['5_1'][which][0], rows['7_1'][which][0]]); e = 0.5*math.hypot(rows['5_1'][which][1], rows['7_1'][which][1])
    l = np.mean([rows['5_1']['law_' + which], rows['7_1']['law_' + which]])
    print('mean q=5,7 %s: slope %+.5f +- %.5f, law %+.5f, ratio %.2f, z %+.1f' % (which, s, e, l, s/l, (s-l)/e))
json.dump(rows, open('slopes.json', 'w'), indent=1)
