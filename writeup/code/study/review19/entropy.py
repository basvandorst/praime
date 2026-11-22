import json, math, numpy as np, sympy
r=json.load(open('/home/bas/praime/writeup/data/study/r64.json'))
d={int(k):v for k,v in json.load(open('full.json'))['data'].items()}; N=sum(d.values())
cnt=np.array(list(d.values()),float); pp=cnt/N
Hd=-np.sum(pp*np.log2(pp))+(len(cnt)-1)/(2*N*math.log(2)); print('H_data MM',Hd,'plug-in',-np.sum(pp*np.log2(pp)))
# C2 and K_tw
ps=np.array(list(sympy.primerange(3,2*10**7)),float)
C2=math.exp(np.sum(np.log1p(-1/(ps-1)**2)))
p5=ps[ps>=5]; f=[1-2/p5,1-3/p5,1-4/p5]; mean=(f[0]+2*f[1]+(p5-3)*f[2])/p5
w=[x/mean for x in f]; K=float(np.sum((w[0]*np.log(w[0])+2*w[1]*np.log(w[1])+(p5-3)*w[2]*np.log(w[2]))/p5))
L0,L1=1e11,1.1e11; LAMS=[math.log(L0),math.log((L0+L1)/2),math.log(L1)]
u=np.linspace(LAMS[0],LAMS[-1],401); wt=np.exp(u)/u**2; lam=float(np.sum(u*wt)/np.sum(wt))
rho=2*C2/lam**2; Has=(math.log(math.e/(6*rho))-K)/math.log(2)
print('C2',C2,'K',K,'lam',lam,'H_asym',Has,'data-asym',Hd-Has)
law=np.array([float(np.sum(np.interp(u,LAMS,r['law'][str(k)])*wt)/np.sum(wt)) for k in range(6,2401,6)])
q=1-6*rho
def H(P): P=P[P>0]; return float(-np.sum(P*np.log2(P)))
def geo(m,n=3000): return m*(1-q)*q**np.arange(n)
for DM in (1200,1800,2400):
    L=law[:DM//6]; t=1-L.sum()
    print('DM',DM,'tail mass',t,'H_law',H(np.concatenate([L,geo(t)])),'data-law',Hd-H(np.concatenate([L,geo(t)])))
# data-only test of the shortcut: replace the data beyond D by a geometric tail with ratio q and the same mass
dv=np.zeros(8000//6+1); 
for k,v in d.items(): dv[k//6-1]=v/N
for D in (1200,2400):
    m=dv[D//6:].sum(); print('data: shortcut beyond',D,'mass',m,'H change',H(np.concatenate([dv[:D//6],geo(m)]))-H(dv))
# structured tail beyond 2400: use the discrete leading-order shape S(d) e^{-rho d} rescaled to the remaining mass
dd=np.arange(2406,30001,6); S=np.ones(len(dd))
for p in sympy.primerange(5,30001):
    nu=np.where(dd%p==0,2,np.where((dd%p==2)|(dd%p==p-2),3,4)); S*=(1-nu/p)/(1-2/p)/(1-1/p)**2
for DM in (1200,2400):
    dd2=np.arange(DM+6,30001,6); S2=np.ones(len(dd2))
    for p in sympy.primerange(5,30001):
        nu=np.where(dd2%p==0,2,np.where((dd2%p==2)|(dd2%p==p-2),3,4)); S2*=(1-nu/p)/(1-2/p)/(1-1/p)**2
    L=law[:DM//6]; t=1-L.sum(); sh=S2*np.exp(-rho*dd2); sh*=t/sh.sum()
    print('DM',DM,'structured tail: H_law',H(np.concatenate([L,sh])),'data-law',Hd-H(np.concatenate([L,sh])))
# data structure check: what does the geometric shortcut cost relative to a structured tail of the same mass beyond 2400
