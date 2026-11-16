import numpy as np
from lib import *
N=10**7; s=sieve(N); P=np.nonzero(s)[0]
# --- difference comb: pairs (p, p+h), p+h<=X, via FFT autocorrelation
X=4*10**6; L=2**23; w=np.zeros(L); w[:X+1]=s[:X+1]; F=np.fft.rfft(w); ac=np.fft.irfft(np.abs(F)**2,L); Hm=420
cnt=np.rint(ac[1:Hm+1]).astype(np.int64)
C2=0.6601618158468696
def S(h):
    if h%2: return 0.0
    r=2*C2; x=h
    while x%2==0: x//=2
    p=3
    while p*p<=x:
        if x%p==0:
            r*=(p-1)/(p-2)
            while x%p==0:x//=p
        p+=2
    if x>1: r*=(x-1)/(x-2)
    return r
t=np.arange(2,X+1,dtype=float); 
pred=np.array([S(h)*np.sum(1/(np.log(t[:len(t)-h])*np.log(t[:len(t)-h]+h))) for h in range(1,Hm+1)])
ev=np.arange(2,Hm+1,2)-1
ratio=cnt[ev]/pred[ev]; print("difference comb: h even 2..%d, count/HL ratio mean %.4f sd %.4f min %.4f max %.4f"%(Hm,ratio.mean(),ratio.std(),ratio.min(),ratio.max()))
print("  h=2:",cnt[1],"pred %.0f"%pred[1]," h=6:",cnt[5],"pred %.0f"%pred[5]," h=30:",cnt[29],"pred %.0f"%pred[29]," h=210:",cnt[209],"pred %.0f"%pred[209])
print("  odd h total pairs:",int(cnt[0::2].sum()),"(only those with p=2)")
print("  count(h)/count(2): h=6 %.3f h=30 %.3f h=210 %.3f ; HL 2, 8/3=2.667, 3.2"%(cnt[5]/cnt[1],cnt[29]/cnt[1],cnt[209]/cnt[1]))
img=np.zeros((360,Hm*3,3))+np.array((14,24,40)); mx=cnt.max()
for h in range(2,Hm+1,2):
    hh=int(cnt[h-1]/mx*340); col=(244,114,182) if h%30==0 else (253,230,138) if h%6==0 else (111,163,214)
    img[360-hh:,(h-1)*3:(h-1)*3+4]=col
    ph=int(pred[h-1]/mx*340); img[360-ph-1:360-ph+1,(h-1)*3-1:(h-1)*3+5]=(255,255,255)
png("comb.png",img)
# --- zero spectrogram of psi
lam=np.zeros(N+1); 
for p in P:
    if p>3200: break
    q=p
    while q<=N: lam[q]=np.log(p); q*=p
lam[P[P>3200]]=np.log(P[P>3200]); psi=np.cumsum(lam)
u0,u1=np.log(50),np.log(N); nu=2**15; u=np.linspace(u0,u1,nu); x=np.exp(u); xi=np.floor(x).astype(int)
f=(psi[xi]-x)/np.sqrt(x)
du=u[1]-u[0]; win=int(4.0/du); hop=win//40; hann=np.hanning(win); gam=np.linspace(5,60,440)
cols=[]
for st in range(0,nu-win,hop):
    seg=(f[st:st+win]-f[st:st+win].mean())*hann; uu=u[st:st+win]
    cols.append(np.abs(np.exp(-1j*np.outer(gam,uu))@seg)*du)
Sg=np.array(cols).T; print("spectrogram shape",Sg.shape,"u range %.2f..%.2f"%(u0,u1))
# global spectrum and peaks
hw=np.hanning(nu); G=np.abs(np.exp(-1j*np.outer(gam,u))@((f-f.mean())*hw))*du
pk=[gam[i] for i in range(1,len(gam)-1) if G[i]>G[i-1] and G[i]>G[i+1] and G[i]>0.25*G.max()]
print("  peaks of the whole-range spectrum at t =",np.round(pk,2))
print("  zeta zeros: 14.13 21.02 25.01 30.42 32.94 37.59 40.92 43.33 48.01 49.77 52.97 56.45 59.35")
img=seq(Sg[::-1],hi=np.percentile(Sg,99.7),stops=HOT); png("spectrogram.png",up(img,2) if img.shape[1]<700 else img); print("  img",img.shape)
np.save("spec_G.npy",np.vstack([gam,G]))
