import numpy as np, math
S=2048; D='/home/bas/scratch/study-review20/'
def Hmm(c, mm=True):
    c=np.asarray(c,float); c=c[c>0]; n=c.sum(); p=c/n
    return -np.sum(p*np.log2(p)) + ((len(c)-1)/(2*n*math.log(2)) if mm else 0)
def marg(keys,c):
    u,inv=np.unique(keys,return_inverse=True); return np.bincount(inv,weights=c)
def cmi(k,c,mm=True):
    g1=k//(S*S); g2=(k//S)%S; g3=k%S
    return Hmm(marg(g1*S+g2,c),mm)+Hmm(marg(g2*S+g3,c),mm)-Hmm(marg(g2,c),mm)-Hmm(c,mm)
d=np.load(D+'pairs_full.npz'); k,c=d['k'],d['c'].astype(float)
I=Hmm(marg(k//S,c))+Hmm(marg(k%S,c))-Hmm(c); I0=Hmm(marg(k//S,c),0)+Hmm(marg(k%S,c),0)-Hmm(c,0)
lam=float(np.mean(np.log(np.linspace(1e12,1e12+1e10,1001))))
print('pairs N',c.sum(),'cells',len(c),'MI MM %.6f plugin %.6f'%(I,I0),'lam',lam,'lam*(I-Iinf)',lam*(I-0.2696508475891928))
rng=np.random.default_rng(1)
def shuffle_table(k,c):
    """random table with the (g1,g2) and (g2,g3) margins fixed: g1 and g3 paired at random within each g2"""
    c=c.astype(np.int64); g1=k//(S*S); g2=(k//S)%S; g3=k%S
    outk=[]; outc=[]
    for b in np.unique(g2):
        m=g2==b
        a1=marg(g1[m],c[m]).astype(np.int64); u1=np.unique(g1[m])
        a3=marg(g3[m],c[m]).astype(np.int64); u3=np.unique(g3[m])
        rem=a3.copy()
        for x,r in zip(u1.tolist(),a1.tolist()):
            dr=rng.multivariate_hypergeometric(rem,r,method='marginals') if rem.sum()>r else rem.copy()
            rem-=dr; nz=dr>0
            outk.append((x*S+b)*S+u3[nz]); outc.append(dr[nz])
    return np.concatenate(outk), np.concatenate(outc).astype(float)
for w in (300000000,1000000000,3000000000,10000000000):
    d=np.load(D+'tri_%d.npz'%w); k,c=d['k'],d['c'].astype(float); N=c.sum()
    J=cmi(k,c); J0=cmi(k,c,False)
    bs=[]
    for r in range(3 if w<10**10 else 2):
        ks,cs=shuffle_table(k,c); bs.append(cmi(ks,cs))
    print('W %.0e triples %d cells %d  J_MM %.6f plugin %.6f  shuffle bias %s  N*bias %.0f'%(w,N,len(c),J,J0,np.round(bs,6),N*np.mean(bs)),flush=True)
