import sys,os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math
from scipy import stats as st
from stats import S2_blocks
L=10**8;H=10**6
tot=[]
for x0,f in [(3*10**11,'pr_3e11.npy'),(6*10**11,'pr_6e11.npy'),(9*10**11,'pr_9e11.npy')]:
    pr=np.load(f)
    chi,col,win,win2,Ph,X=S2_blocks(pr,x0,L,H,3,2)
    nw=np.searchsorted(pr,x0+L)
    same=(pr[:nw]%3)==(pr[2:nw+2]%3)
    blk=(pr[:nw]-x0)//H
    Ek=np.bincount(blk,weights=same,minlength=100); N=np.bincount(blk,minlength=100)
    pi=Ek.sum()/N.sum()
    disp=(Ek-pi*N)**2/(N*pi*(1-pi))
    print(x0,'corr(within, E_k dispersion) %.4f'%np.corrcoef(win,disp)[0,1],' sum within/99 %.4f  sum disp/99 %.4f'%(win.sum()/99,disp.sum()/99))
    tot.append(win.sum())
S=sum(tot); dof=297
law=0.993
print('pooled within %.4f'%(S/dof))
print('P(chi2_297 * law/297 <= observed) = %.2e'%st.chi2.cdf(S/law,dof), ' z(normal, law sd) = %.2f'%((S/dof-law)/(law*math.sqrt(2/dof))))
print('one-sided p -> normal equivalent z = %.2f'%st.norm.ppf(st.chi2.cdf(S/law,dof)))
