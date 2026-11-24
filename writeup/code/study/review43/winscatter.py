import numpy as np, json, math, os
hl=json.load(open('hlpred.json'))
W=['2e11','2.5e11','3e11','3.5e11','4e11','4.5e11','5.5e11','6e11','7e11','7.5e11','8e11','9e11','1e12','1.1e12','1.2e12']
D={w:np.load('an_%s.npz'%w) for w in W}
L=1e8
print('S1 over 15 windows 2e11-1.2e12 (12 fresh + 3 registered); window-scatter chi2 = sum((m_w - law_w)/e_w)^2 on 15')
for q in [3,4,5,7,8,9,12,30,210]:
    for H in [10**4,10**5,10**6,10**7]:
        ms=[];es=[];ls=[]; allres=[]
        for w in W:
            v=D[w]['S1_%d_%d'%(q,H)]; lam=math.log(float(w)+L/2)
            law=1+hl['%d_%d'%(q,H)]['a_exact']/lam
            ms.append(v.mean()); es.append(v.std(ddof=1)/math.sqrt(len(v))); ls.append(law)
            allres.append(v-law)
        ms=np.array(ms);es=np.array(es);ls=np.array(ls)
        allres=np.concatenate(allres)
        m=np.mean(ms); l=np.mean(ls)
        err_blocks=allres.std(ddof=1)/math.sqrt(len(allres))
        err_windows=np.std(ms-ls,ddof=1)/math.sqrt(len(W))
        chi=np.sum(((ms-ls)/es)**2)
        print('q=%3d H=%.0e meas %.4f law %.4f  diff %+.4f  err(blocks) %.4f err(windows) %.4f  z %+.2f  window chi2 %.1f/15'%(q,H,m,l,m-l,err_blocks,err_windows,(m-l)/max(err_blocks,err_windows),chi))
