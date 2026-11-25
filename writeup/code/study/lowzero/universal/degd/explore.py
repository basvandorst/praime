import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math, time
from degd_lib import *
cases = {
 'X5':  (5.0, [1.0], [(400,440,4.0)]),
 'D12': (1.0, [5.5,6.5], [(200,260,5.0)]),
 'S11a':(121.0,[1.0,1.0,2.0], [(40,50,4.0),(60,80,5.0),(30,36,3.0),(36,44,4.0)]),
}
for key,(Q,mu,wins) in cases.items():
    for (T1,T2,D) in wins:
        t0=time.time()
        tg,h = grid(T1,T2,D,h=0.05)
        f = H_box(tg,T1,T2,D)*np.exp(2j*theta(tg,Q,mu))
        lu = np.linspace(0.5, 22, 4000)
        I = osc_sum(lu, tg, f, h)
        w = np.abs(I)*np.exp(lu/2)
        big = lu[w>1e-11]
        print(key,(T1,T2,D),'ln u with |I|sqrt(u)>1e-11:', round(big.min(),2), round(big.max(),2), 'u_max=%.3g'%math.exp(big.max()), 'K=%d'%len(tg), 'secs %.1f'%(time.time()-t0))
        print('   2theta\'(ends, T1, T2)=', np.round(2*theta1(np.array([tg[0],T1,T2,tg[-1]]),Q,mu),2))
