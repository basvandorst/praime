import sys; sys.path.insert(0,'/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/referee_al')
from s3_level import Dq, DN, nup, pari
import numpy as np, math
from fractions import Fraction as Fr
print('y     Dq      DN     ++      -+')
for y in (0.05,0.1,0.2,0.5,1,1.5,2,3):
    a,b=Dq(y),DN(y); print('%4.2f %7.3f %7.3f %7.3f %7.3f'%(y,a,b,a+b,a-b))
# Martin figure: M^eps(x;2) = int_1^2 D^eps(x/t) t dt / 1.5, X = 6000
ts=np.linspace(1,2,2001)
for ell in (250,800,2000,4000):
    x=ell/6000; out=[]
    for e1,eq in ((1,1),(-1,-1),(1,-1),(-1,1)):
        f=np.array([(eq*Dq(x/t)+e1*eq*DN(x/t))*t for t in ts]); out.append(np.trapezoid(f,ts)/1.5/4)
    print('l=%d M/4 (++,--,+-,-+):'%ell, ' '.join('%6.3f'%v for v in out))
# eigenspace dims at N = 20206 from Martin's formula at l = 1 (s=0 only for Q = q, N; Q = 2: s = 0, +-1)
q=10103; N=2*q
H=lambda n: Fr(str(pari.qfbhclassno(n)))
def xi(D,l):
    return l*H(-D//(l*l))/H(-D)-1 if D%(l*l)==0 else int(pari.kronecker(D,l))-1
trq = -xi(-4*q,2)*H(4*q)/2 - 1
trN = -H(4*N)/2 + 1
tr2 = -sum(xi(s*s*4-8, q)*H(8-4*s*s) for s in (-1,0,1))/2 - 1
dim=int(pari('mfdim([%d,2],0)'%N))
print('N=20206 dim',dim,'trW_2',tr2,'trW_q',trq,'trW_N',trN, 'eig dims', [ (dim+e1*tr2+eq*trq+e1*eq*trN)/4 for e1 in (1,-1) for eq in (1,-1)], 'sqrt(N)=%.0f'%math.sqrt(N))
