import numpy as np, math, sys
K=int(sys.argv[1]) if len(sys.argv)>1 else 21
N=1<<K
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
s=sieve(N)[:N]
def fwht(a):
    a=a.copy(); h=1; n=len(a)
    while h<n:
        a=a.reshape(-1,2,h); x=a[:,0,:]+a[:,1,:]; y=a[:,0,:]-a[:,1,:]; a=np.stack([x,y],axis=1).reshape(-1); h*=2
    return a
def spectrum(ind):
    f=ind[1::2].astype(float); mu=f.mean(); W=fwht(f-mu)/len(f); return W/mu,mu
W,mu=spectrum(s)
print("odd numbers below 2^%d: share of primes %.4f. For a set of bit layers (a mask), split the odd numbers by the parity of their bits in those layers;"%(K,mu))
print("bias = (share of primes in the even half - share in the odd half) / (2 x overall share). Noise for one mask: about %.4f"%(math.sqrt((1-mu)/mu/(N/2))))
o=np.argsort(-np.abs(W))[:28]
def pat(m): return format(int(m)<<1,'0%db'%K)
print("top masks (layers written from the highest bit to bit 0; bit 0 is always 1 and left out):")
for m in o: print("   %s  layers %2d  bias %+.4f"%(pat(m),bin(int(m)).count('1'),W[m]))
ab=np.sort(np.abs(W))[::-1]; e=(W**2).sum()
print("energy: sum of squared biases %.4f (= (1-share)/share for a set with no structure: %.4f / %d masks x %d = %.4f); top 10 masks hold %.3f of it, top 100 %.3f, top 1000 %.3f"%(e,(1-mu)/mu,len(W),len(W),(1-mu)/mu,(ab[:10]**2).sum()/e,(ab[:100]**2).sum()/e,(ab[:1000]**2).sum()/e))
np.save('W%d.npy'%K,W)
def model(ps):
    ind=np.ones(N,bool); ind[:2]=False
    for p in ps: ind[::p]=False
    return ind
for ps in ((3,),(5,),(7,),(3,5),(3,5,7),(3,5,7,11,13)):
    Wm,_=spectrum(model(ps)); print("numbers with no factor %s: correlation of its biases with those of the primes over all masks %.3f; over the primes' top 1000 masks %.3f"%(ps,np.corrcoef(W,Wm)[0,1],np.corrcoef(W[np.argsort(-np.abs(W))[:1000]],Wm[np.argsort(-np.abs(W))[:1000]])[0,1]))
