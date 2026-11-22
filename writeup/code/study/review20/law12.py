import numpy as np, math, sys
sys.path.insert(0,'/home/bas/scratch/study-review20'); from sv import base_primes
lam=float(np.mean(np.log(np.linspace(1e12,1e12+1e10,1001))))
def H(w):
    w=w[w>0]; p=w/w.sum(); return -np.sum(p*np.log2(p))
def margH(keys,w):
    u,inv=np.unique(keys,return_inverse=True); return H(np.bincount(inv,weights=w))
def pair_law(G):
    ev=np.arange(2,G+1,2); A,B=np.meshgrid(ev,ev,indexing='ij'); A=A.ravel(); B=B.ravel(); m=A+B<=G; A,B=A[m],B[m]
    lw=-(A+B)/lam; alive=np.ones(len(A),bool)
    for p in base_primes(G)[1:].tolist():
        r1=A%p; r2=(A+B)%p
        nu=1+(r1!=0)+((r2!=0)&(r2!=r1))
        f=1-nu/p; alive&=f>0; lw+=np.log(np.where(f>0,f,1))
    w=np.where(alive,np.exp(lw-lw.max()),0)
    return margH(A,w)+margH(B,w)-H(w)
def tri_law(G):
    ev=np.arange(2,G+1,2,dtype=np.int32)
    A,B=np.meshgrid(ev,ev,indexing='ij'); A=A.ravel(); B=B.ravel(); m=A+B<=G-2; A,B=A[m],B[m]
    As=[];Bs=[];Cs=[]
    for a,b in zip(A.tolist(),B.tolist()):
        c=np.arange(2,G-a-b+1,2,dtype=np.int32); As.append(np.full(len(c),a,np.int32)); Bs.append(np.full(len(c),b,np.int32)); Cs.append(c)
    A=np.concatenate(As); B=np.concatenate(Bs); C=np.concatenate(Cs)
    lw=-(A+B+C).astype(float)/lam; alive=np.ones(len(A),bool)
    for p in base_primes(G)[1:].tolist():
        r1=A%p; r2=(A+B)%p; r3=(A+B+C)%p
        nu=1+(r1!=0)+((r2!=0)&(r2!=r1))+((r3!=0)&(r3!=r1)&(r3!=r2))
        f=1-nu/p; alive&=f>0; lw+=np.log(np.where(f>0,f,1))
    w=np.where(alive,np.exp(lw-lw.max()),0); S=1024
    k12=A*S+B; k23=B*S+C
    return margH(k12,w)+margH(k23,w)-margH(B,w)-H(w), len(A)
print('lam',lam)
for f in (30,40):
    G=int(f*lam)//2*2; print('pair law G=%d (%d lam): I = %.6f'%(G,f,pair_law(G)),flush=True)
for f in (20,):
    G=int(f*lam)//2*2; J,n=tri_law(G); print('triple law G=%d: J = %.6f (%d triples)'%(G,J,n),flush=True)
