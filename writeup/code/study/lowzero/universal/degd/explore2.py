import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math, time
from degd_lib import *
cases = [('X5',5.0,[1.0],(400,440,4.0)),('D12',1.0,[5.5,6.5],(200,260,5.0)),('S11a',121.0,[1.0,1.0,2.0],(40,50,4.0)),('S11a',121.0,[1.0,1.0,2.0],(36,44,4.0))]
for key,Q,mu,(T1,T2,D) in cases:
    tg,h = grid(T1,T2,D,h=0.05)
    f = H_box(tg,T1,T2,D)*np.exp(2j*theta(tg,Q,mu))
    lu = np.arange(1,22,0.5)
    I = osc_sum(lu, tg, f, h)
    print(key,(T1,T2,D))
    print('  '+' '.join('%.1f:%.0e'%(a,b) for a,b in zip(lu,np.abs(I))))
