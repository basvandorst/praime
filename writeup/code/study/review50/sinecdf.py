import numpy as np, json
from numpy.polynomial.legendre import leggauss
def E0(s,m=40):
    x,w=leggauss(m); x=(x+1)*s/2; w=w*s/2
    d=x[:,None]-x[None,:]
    K=np.sinc(d)  # sin(pi d)/(pi d)
    A=np.sqrt(w)[:,None]*K*np.sqrt(w)[None,:]
    return np.linalg.det(np.eye(m)-A)
def cdf(s0,h=1e-5):
    # F(s0)=1+E'(s0)
    return 1+(E0(s0+h)-E0(s0-h))/(2*h)
st=json.load(open('base_stats.json'))
tot={k:0 for k in ['n','n05','n10','n15','n20','n30']}
for name in ['B','Dn','Fn','C','Bp']:
    for x in st[name]:
        for k in tot: tot[k]+=x[k]
print('total spacings',tot['n'])
for s0,k in [(0.05,'n05'),(0.1,'n10'),(0.15,'n15'),(0.2,'n20'),(0.3,'n30')]:
    F=cdf(s0); exp=F*tot['n']
    print('s<%.2f  observed %6d  sine expected %8.1f  ratio %.3f +- %.3f   (pi^2/9 s^3 approx %.2e vs exact %.2e)'%(s0,tot[k],exp,tot[k]/exp,np.sqrt(tot[k])/exp,np.pi**2/9*s0**3,F))
# per set ratio at 0.3 and at 0.1
for name in ['B','Dn','Fn','C','Bp','R120']:
    n=sum(x['n'] for x in st[name]); 
    print(name, ' '.join('%s %.3f'%(k,sum(x[k] for x in st[name])/(cdf(s0)*n)) for s0,k in [(0.05,'n05'),(0.1,'n10'),(0.15,'n15'),(0.2,'n20'),(0.3,'n30')]))
