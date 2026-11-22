import json, math, numpy as np
r=json.load(open('/home/bas/praime/writeup/data/study/r64.json'))
d={int(k):v for k,v in json.load(open('full.json'))['data'].items()}; N=sum(d.values())
L0,L1=1e11,1.1e11; LAMS=[math.log(L0),math.log((L0+L1)/2),math.log(L1)]
u=np.linspace(LAMS[0],LAMS[-1],401); wt=np.exp(u)/u**2
law={int(k):float(np.sum(np.interp(u,LAMS,v)*wt)/np.sum(wt)) for k,v in r['law'].items()}
lam=float(np.sum(u*wt)/np.sum(wt)); print('mean lambda',lam, r['compare']['mean_lambda'])
bands=[(lo,lo+24) for lo in range(6,600,30)]+[(lo,lo+54) for lo in range(606,1200,60)]
chi=0; zs=[]
for lo,hi in bands:
    c=sum(v for k,v in d.items() if lo<=k<=hi); e=N*sum(law.get(k,0) for k in range(lo,hi+1,6))
    chi+=(c-e)**2/e; zs.append((c-e)/math.sqrt(e))
print('chi2',chi,'max|z|',max(map(abs,zs)), 'stored',r['compare']['chi2'])
m1200=sum(v for k,v in law.items() if k<=1200); mall=sum(law.values())
tc=sum(v for k,v in d.items() if k>1200)
te=N*(1-m1200); print('mass<=1200',m1200,'stored',r['compare']['law_mass_to_1200'],'tail',tc,te,'rel',tc/te-1,'z',(tc-te)/math.sqrt(te))
# better variance: binomial
print('z binomial', (tc-te)/math.sqrt(N*(1-m1200)*m1200))
print('mass<=2400',mall,'rerun tail_law would be',N*(1-mall),'vs tail_count',tc)
t2=sum(v for k,v in d.items() if k>2400); print('data >2400',t2, 'law >2400', N*(1-mall), 'z',(t2-N*(1-mall))/math.sqrt(N*(1-mall)))
# 1206-2400 band
c=sum(v for k,v in d.items() if 1200<k<=2400); e=N*(mall-m1200); print('1206-2400',c,e,c/e-1,(c-e)/math.sqrt(e))
# finer: bands of 120 in 1206..2400
for lo in range(1206,2400,240):
    c=sum(v for k,v in d.items() if lo<=k<lo+240); e=N*sum(law.get(k,0) for k in range(lo,lo+240,6)); print(lo,c,round(e),round((c-e)/math.sqrt(e),2))
# exponential baseline check
pi6=6*1.3203236316937/lam**2
chiP=sum(((sum(v for k,v in d.items() if lo<=k<=hi))-N*sum(pi6*(1-pi6)**(k//6-1) for k in range(lo,hi+1,6)))**2/(N*sum(pi6*(1-pi6)**(k//6-1) for k in range(lo,hi+1,6))) for lo,hi in bands)
print('chiP',chiP)
