import numpy as np, math
K=24; N=1<<K; n=np.arange(N,dtype=np.int64)
def pmodarr(f):
    r=n.copy(); df=f.bit_length()-1
    for b in range(K-1,df-1,-1):
        sel=((r>>b)&1)==1; r[sel]^=f<<(b-df)
    return r
print("all numbers below 2^L: how often do 'M = 2^w - 1 divides n' (with carries) and 'the strip of w ones divides n' (without carries) hold together, against chance")
print("  w   M     " + " ".join("L=%-6d"%L for L in (12,16,20,24)))
res={}
for w in (2,3,4,5,6,7,8):
    M=(1<<w)-1; x=(n%M==0); y=(pmodarr(M)==0); row=[]
    for L in (12,16,20,24):
        T=1<<L; a=int(x[:T].sum()); b=int(y[:T].sum()); c=int((x[:T]&y[:T]).sum()); row.append(c/(a*b/T))
    res[w]=row; print("  %d  %-5d "%(w,M)+" ".join("%.4f  "%v for v in row))
print("controls (different widths): ")
for (p,f) in ((3,7),(7,3),(7,31),(31,7),(5,7),(9,7),(7,5)):
    x=(n%p==0); y=(pmodarr(f)==0); T=N; print("  %d | n with carries, pattern %d | n without: ratio %.4f"%(p,f,int((x&y).sum())/(int(x.sum())*int(y.sum())/T)))
