# k = 1 with a smooth weight: sum_gamma w(gamma)|L'(rho)|^2 against int w(t) m1(t) dt (m1 = dS/dT of the Gonek main term)
import numpy as np, json, sys
sys.path.insert(0,'.')
from gonek_k1 import poly_coeffs
import sympy
z=np.load('../all.npz'); d=json.load(open('../perL.json'))
def w(t,a=100,b=600,tap=50):
    t=np.asarray(t,float); out=np.ones_like(t)
    m=t<a+tap; out[m]=0.5-0.5*np.cos(np.pi*(t[m]-a)/tap)
    m=t>b-tap; out[m]=0.5-0.5*np.cos(np.pi*(b-t[m])/tap)
    out[(t<a)|(t>b)]=0; return out
xt,wt=np.polynomial.legendre.leggauss(200); tt=(xt+1)/2*500+100; wtt=wt/2*500
res=[]
for r in d:
    q=r['q']; P=[float(x) for x in poly_coeffs(q)]
    x=np.log(q*tt/(2*np.pi)); m1=sum(P[m]*(x**m+m*x**(m-1 if m>0 else 0)*(m>0)) for m in range(5))/(2*np.pi)
    pred_w=np.sum(wtt*w(tt)*m1); pred_s=np.sum(wtt*m1)
    g=z['%s_%d_g'%(r['set'],r['D'])]; y=z['%s_%d_y'%(r['set'],r['D'])]
    res.append((r['D'],np.log(np.sum(w(g)*np.exp(2*y))/pred_w), np.log(r['S1']/pred_s), np.log(r['S1']/r['S1pred'])))
a=np.array(res); keep=np.abs(a[:,0])!=1153
print('check sharp integral vs S(600)-S(100): max |diff| in log', np.abs(a[:,2]-a[:,3]).max())
for lab,c in [('sharp',2),('smooth taper 50',1)]:
    v=a[keep,c]; print('%-16s mean log ratio %+.5f +- %.5f ; sd per L-function %.5f ; max |.| %.4f'%(lab,v.mean(),v.std(ddof=1)/np.sqrt(len(v)),v.std(ddof=1),np.abs(v).max()))
