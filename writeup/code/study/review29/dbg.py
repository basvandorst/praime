from rz import *
g=load(2419); d=0.1; c=coef(g,d)
for span in (14.0,):
  for h in (1e-4,5e-5):
    L1=math.log(1e8); u=np.arange(L1-span,L1+h/2,h); u[-1]=L1
    S=S_grid(g,c,u); f=np.exp(u-2*L1)*S**2
    trap=h*(f.sum()-0.5*f[0]-0.5*f[-1]); 
    print(h,len(u),'trap',2*trap/d, 'last step',u[-1]-u[-2])
print('direct',R_direct(g,d,1e8))
