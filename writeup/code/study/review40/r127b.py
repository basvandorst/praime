import numpy as np
exec(open('r127.py').read().split("t0 = time.time()")[0])
counts = np.load('counts127.npy')
g = np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_40000_flint.txt'); rho = 0.5+1j*g
def lisum(x):
    z = rho*np.log(x); return np.sum(np.exp(z)/z*(1+1/z+2/z**2+6/z**3+24/z**4))
m1=[];m2=[];z1=[];sh=[]
for j in range(NB):
    a = X0+j*BLK; xs = a+H*np.arange(BLK//H+1)
    S = np.array([lisum(float(x)) for x in xs]); dev = -2*np.real(S[1:]-S[:-1])
    main = H/np.log(xs[:-1]+H/2)
    r = counts[j]-main
    m1.append(np.var(r, ddof=1)/main.mean()); m2.append(np.mean(r**2)/main.mean())
    xm=a+BLK/2; lx=np.log(xm); h=np.log1p(H/xm)
    gg=np.linspace(g[-1],2e7,4000001)
    tail=np.trapezoid(8*xm*np.sin(gg*h/2)**2/(gg**2*lx**2)*np.log(gg/(2*np.pi))/(2*np.pi),gg)+4*xm/lx**2*np.log(2e7/(2*np.pi))/(2*np.pi)/2e7
    z1.append((np.var(dev,ddof=1)+tail)/main.mean()); sh.append(tail/(np.var(dev,ddof=1)+tail))
print('measured var(counts-H/ln mid) ddof1 %.4f, mean square %.4f; zeros %.4f tail share %.3f'%(np.mean(m1),np.mean(m2),np.mean(z1),np.mean(sh)))
print('block sd of measured %.4f'%np.std(m1,ddof=1))
