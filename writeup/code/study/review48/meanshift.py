# ratios-conjecture pair correlation -> per-L-function mean shift of log|L'(rho)| from the conductor's primes,
# diagonal (first order) + off-diagonal (swapped term) exactly, relative to the q = 1 arithmetic.
import numpy as np, json, sys
g=np.load(sys.argv[1]); r=g['r']; w=g['w']; Z=g['Z']
R0s=[15.,30.,50.]
d=json.load(open('../perL.json'))
xt,wt=np.polynomial.legendre.leggauss(24); tt=(xt+1)/2*500+100; wtt=wt/2*500
lnr=np.log(r)
def rho_minus1(bad):
    s=np.zeros(len(r),complex)
    for p in bad:
        lp=np.log(p); dl=(4/p)*np.sin(r*lp/2)**2/((1-2/p)+np.exp(-1j*r*lp)/p)
        s+=np.log1p(dl)
    return np.expm1(s)
def off(bad,Ls,R0):
    h=Z*rho_minus1(bad)*np.exp(-(r/R0)**2)
    # I(L) = 2 int_0^inf ln r Re[e^{-irL} h(r)] dr
    return np.array([2*np.sum(w*lnr*np.real(np.exp(-1j*r*L)*h)) for L in Ls])
# single-prime checks against the asymptotic -h0/L^2 (h0 = ln^2 p/(p-1))
for p in [2,3,5,7,13]:
    for L in [8.,10.,12.]:
        I=[off([p],[L],R0)[0] for R0 in R0s]
        first=np.log(p)*np.log(p/(p-1))/L
        print('p=%2d L=%4.1f first %.5f  off(exact,R0=15/30/50) %s  asym -h0/L^2 %.5f  ratio total/first %.4f'%(
            p,L,first,' '.join('%.5f'%(x/(np.pi*L)) for x in I),-np.log(p)**2/(p-1)/L**2,(first+I[-1]/(np.pi*L))/first))
out=[]
for rr in d:
    q=rr['q']; Ls=np.log(q*tt/(2*np.pi)); wz=wtt*Ls; wz/=wz.sum()
    first=sum(np.log(p)*np.log(p/(p-1)) for p in rr['bad'])*np.sum(wz/Ls)
    offs=[np.sum(wz*off(rr['bad'],Ls,R0)/(np.pi*Ls)) for R0 in R0s]
    out.append(dict(D=rr['D'],first=first,off=offs))
json.dump(out,open(sys.argv[2],'w'))
print('done',len(out))
