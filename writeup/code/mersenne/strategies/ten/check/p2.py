import numpy as np
from lib import *
N=10**7; s=sieve(N); P=np.nonzero(s)[0]
# --- gap return map
Q=P[P>=5]; g=np.diff(Q); G=40
a=g[:-1]//2; b=g[1:]//2; m=(a<=G)&(b<=G)
H=np.zeros((G+1,G+1)); np.add.at(H,(a[m],b[m]),1)
ma=np.bincount(g//2,minlength=G+1)[:G+1].astype(float); E=np.outer(ma,ma)/len(g)
forb=[(i,j) for i in range(1,G+1) for j in range(1,G+1) if (2*i)%6==(2*j)%6 and (2*i)%6!=0]
print("return map: forbidden cells (both gaps = 2 or both = 4 mod 6):",len(forb),"of",G*G,"; pairs in them:",int(sum(H[i,j] for i,j in forb)))
R=np.where(E>20,H/np.maximum(E,1e-9),np.nan)
d=np.array([R[i,i] for i in range(3,G+1,3) if not np.isnan(R[i,i])]); print("diagonal (g,g) with 6|g: obs/indep ratios",np.round(d[:8],3))
allowed=np.array([[not((2*i)%6==(2*j)%6 and (2*i)%6!=0) for j in range(G+1)] for i in range(G+1)])
print("mean ratio allowed cells (E>20): %.3f ; min %.3f max %.3f"%(np.nanmean(R[allowed&~np.isnan(R)]),np.nanmin(R[allowed&~np.isnan(R)&(R>0)]),np.nanmax(R)))
img=div(np.nan_to_num(np.log2(np.maximum(R,1/8)),nan=0.0),3.0); img[np.isnan(R)]=(14,24,40)
png("returnmap.png",up(img[1:,1:][::-1],16))
# --- diffraction
L=30030*64; M=8*10**6; w=s[M:M+L].astype(float); n=w.sum(); F=np.abs(np.fft.rfft(w))**2/n
print("diffraction: window",M,"+",L,"primes",int(n))
for q,aa in [(2,1),(3,1),(5,1),(5,2),(6,1),(7,1),(10,1),(15,1),(30,1),(4,1),(9,1),(8,1),(11,1),(210,1)]:
    from math import gcd
    k=aa*L/q; k0=int(round(k)); pk=F[k0]
    ph=q; x=q
    for p in (2,3,5,7,11,13):
        if x%p==0:
            ph=ph//p*(p-1)
            while x%p==0:x//=p
    sqf=all(q%(p*p) for p in (2,3,5,7))
    print("  a/q=%d/%d  S/N=%.5f   mu^2/phi^2=%.5f"%(aa,q,pk/n,(1/ph**2 if sqf else 0)))
print("  median background S: %.3f (Poisson would be ~1-density)"%np.median(F[1:]))
# picture: S(k) on log scale, k/L in [0,1/2], max-pooled to 1400 columns
Fh=F[1:L//2+1]; cols=1400; pool=Fh[:(len(Fh)//cols)*cols].reshape(cols,-1).max(axis=1)
y=np.log10(pool); print("  log10 S range",y.min(),y.max())
img=np.zeros((420,cols,3))+np.array((14,24,40)); yy=((y-(-1))/(np.log10(n)+0.2+1)*419).astype(int).clip(0,419)
for c in range(cols): img[419-yy[c]:,c]=(111,163,214) if yy[c]<140 else (244,114,182)
png("diffraction.png",img)
# --- 2D diffraction of rows layout width 1024
W=1024; A=s[1:W*W+1].reshape(W,W).astype(float) if N>=W*W else None
if A is not None:
    F2=np.fft.fftshift(np.abs(np.fft.fft2(A-A.mean()))**2); png("diffraction2d_rows.png",seq(np.log10(F2+1),stops=HOT))
# ulam spiral 1001x1001
def ulam(side):
    A=np.zeros((side,side),dtype=np.int64); x=y=side//2; A[y,x]=1; nn=1; step=1
    while nn<side*side:
        for dx,dy in ((1,0),(0,-1)):
            for _ in range(step):
                if nn>=side*side:break
                x+=dx;y+=dy;nn+=1;A[y,x]=nn
        step+=1
        for dx,dy in ((-1,0),(0,1)):
            for _ in range(step):
                if nn>=side*side:break
                x+=dx;y+=dy;nn+=1;A[y,x]=nn
        step+=1
    return A
U=ulam(1023); np.save("ulam1023.npy",U); A=s[U].astype(float); F2=np.fft.fftshift(np.abs(np.fft.fft2(A-A.mean()))**2)
png("diffraction2d_ulam.png",seq(np.log10(F2+1),lo=2.5,stops=HOT)); print("ulam fft log range",np.log10(F2+1).min(),np.log10(F2+1).max())
