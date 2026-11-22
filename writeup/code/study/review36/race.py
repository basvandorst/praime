import json, numpy as np, sympy
from scipy.stats import pearsonr, spearmanr
d=json.load(open('data.json'))
D=np.array([o['D'] for o in d])
P=np.array(list(sympy.primerange(2,100001)))
# chi table: rows D, cols p
chi=np.zeros((len(D),len(P)),dtype=np.int8)
m8=D%8
chi[:,0]=np.where(D%2==0,0,np.where((m8==1)|(m8==7),1,-1))
for j,p in enumerate(P[1:],1):
    qr=np.full(p,-1,np.int8); qr[(np.arange(1,p)**2)%p]=1; qr[0]=0
    chi[:,j]=qr[D%p]
cs=np.cumsum(chi,axis=1,dtype=np.int32)   # N+ - N-
valid=chi!=0
actual=np.array([( (cs[i]<0)&valid[i]).sum()/valid[i].sum() for i in range(len(D))])
lnP=np.log(P.astype(float))
pred=np.zeros(len(D)); pred1=np.zeros(len(D))
for i,o in enumerate(d):
    z=np.array([x for x in o['z'] if 0<x<=30]); rho=0.5+1j*z
    ph=np.exp(1j*np.outer(lnP,z))/rho
    E=1+2*ph.real.sum(axis=1); E1=1+2*ph[:,0].real
    pred[i]=((E>0)&valid[i]).sum()/valid[i].sum(); pred1[i]=((E1>0)&valid[i]).sum()/valid[i].sum()
np.save('race.npy',np.vstack([D,actual,pred,pred1]))
print('Pearson all',pearsonr(pred,actual)[0],'mean abs err',np.abs(pred-actual).mean(),'lowest',pearsonr(pred1,actual)[0], np.abs(pred1-actual).mean())
for s,n in ((D<0,'neg'),(D>0,'pos')):
    print(n, pearsonr(pred[s],actual[s])[0], pearsonr(pred1[s],actual[s])[0])
print('mean share', actual.mean())
for x in (-163,-4,-3,5,-2787):
    i=list(D).index(x); print(x, actual[i],pred[i],pred1[i], d[i]['z'][0])
