# Closed-form column-part prediction: col/dof = tr(pinv(Mt) (I-P')K(I-P))/(C-1), K = HL class covariance / mean.
import sys,os,json,math,numpy as np
from sympy import factorint, totient
hl=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'hlpred.json')))
L=10**8
wins=[float(s) for s in sys.argv[1].split(',')]
def Kmat(q,H,lam):
    A=np.array(hl['%d_%d'%(q,H)]['A']); return np.eye(len(A))+A/lam
def colpred(P,K):
    C=len(P); Mt=sum(np.diag(P[i])-np.outer(P[i],P[i]) for i in range(C))
    I=np.eye(C)
    V=(I-P.T)@K@(I-P)
    return np.trace(np.linalg.pinv(Mt)@V)/(C-1)
for q,k in [(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(10,1),(10,2),(10,3),(5,1),(7,1),(7,2)]:
    for H in [10**4,10**5,10**6]:
        meas=[];pred=[];allv=[]
        for x0 in wins:
            f='an_%s.npz'%('%g'%x0).replace('e+','e')
            d=np.load(f)
            key='S2col_%d_%d_%d'%(q,k,H)
            if key not in d.files: break
            v=d[key]; nb=len(v); allv.append(v*nb/(nb-1)); meas.append(v.mean()*nb/(nb-1))
            P=d['Ph_%d_%d'%(q,k)]
            lam=math.log(x0+L/2)
            pred.append(colpred(P,Kmat(q,H,lam)))
        if not meas: continue
        allv=np.concatenate(allv)
        err=allv.std(ddof=1)/math.sqrt(len(allv))
        print('(%2d,%d) H=%.0e  col meas %.4f +- %.4f   closed form %.4f   z %+.2f'%(q,k,H,np.mean(meas),err,np.mean(pred),(np.mean(meas)-np.mean(pred))/err))
d=np.load('an_3e11.npz')
for k in [1,2,3,4]:
    P=d['Ph_3_%d'%k]; print('q=3 k=%d P='%k, P.round(4).tolist())
