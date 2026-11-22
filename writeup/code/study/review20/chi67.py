import json, math, numpy as np
from scipy.stats import chi2
d=json.load(open('/home/bas/praime/writeup/data/study/r67.json')); cnt={int(k):v for k,v in json.load(open('/home/bas/scratch/study-review20/tw13_counts.json')).items()}
L0,L1=1e13,1e13+1e10; lams=[math.log(L0),math.log((L0+L1)/2),math.log(L1)]
# weight lam over the window by twin density 1/ln^2 t (dt = e^u du)
u=np.linspace(lams[0],lams[-1],2001); wt=np.exp(u)/u**2
law={int(k):float(np.sum(np.interp(u,lams,v)*wt)/np.sum(wt)) for k,v in d['law'].items()}
N=sum(cnt.values()); mass=sum(law.values())
bands=[(lo,lo+24) for lo in range(6,600,30)]+[(lo,lo+54) for lo in range(606,1200,60)]
ch=0; zs=[]; resid=0
for lo,hi in bands:
    c=sum(v for k,v in cnt.items() if lo<=k<=hi); e=N*sum(law.get(k,0) for k in range(lo,hi+1,6)); ch+=(c-e)**2/e; zs.append((c-e)/math.sqrt(e)); resid+=c-e
tc=sum(v for k,v in cnt.items() if k>1200); te=N*(1-mass)
print('N',N,'chi2 %.2f p(chi2_30>=)=%.3f max|z| %.2f'%(ch,chi2.sf(ch,30),max(map(abs,zs))))
print('tail data %d law %.1f diff %.1f ; -(sum band residuals) %.1f'%(tc,te,tc-te,-resid))
q=1-mass; print('tail share %.4f, z Poisson %.3f, z binomial %.3f'%(tc/N,(tc-te)/math.sqrt(te),(tc-te)/math.sqrt(N*q*(1-q))))
lm=float(np.sum(u*wt)/np.sum(wt)); print('mean lam',lm)
# heuristic truncation scaling: overestimate of the law mass below 1200 ~ gamma(4, 1200 rho_tw), rho_tw = 2 C2 / lam^2
from scipy.special import gammainc
def G4(x): return gammainc(4,x)*6
C2x2=1.3203236316937
f11=4001/20498567; l11=25.3768; l13=lm; l9=21.07
for name,l in (('1e9',l9),('1e11',l11),('1e13',l13)): print(name,'1200*rho',1200*C2x2/l**2,'G4',G4(1200*C2x2/l**2))
pred=f11*G4(1200*C2x2/l13**2)/G4(1200*C2x2/l11**2)*N
print('predicted truncation excess at 1e13 if unchanged in scaled units: %.0f gaps, z %.2f'%(pred,pred/math.sqrt(N*q*(1-q))))
print('1e9 check: predicted fraction from 1e11 %.2e vs observed %.2e'%(f11*G4(1200*C2x2/l9**2)/G4(1200*C2x2/l11**2), 989/2963534))
