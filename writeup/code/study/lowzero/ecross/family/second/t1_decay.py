import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, time, numpy as np
from smlib import *
N, k = 10007, 4
for T in [6.0, 12.0, 20.0]:
    t0 = time.time()
    w1 = Window(N, k, T, h=0.02); w2 = Window(N, k, T, h=0.05)
    ts = np.linspace(0.3*T, 3.5*T, 60)
    u = N*(ts/(2*math.pi))**2
    I1 = w1.I(u); I2 = w2.I(u)
    amp = np.abs(-(1/(2*math.pi))*u**-0.5*I1*np.log(u))
    print('T=%g intH=%.4f L=%.4f intH2=%.4f; max|I(h=.02)-I(h=.05)| = %.1e  (%.1fs)' % (T, w1.intH, w1.L, w1.intH2, np.max(np.abs(I1-I2)), time.time()-t0))
    for t, a, i in zip(ts[::3], amp[::3], I1[::3]): print('   t_n=%6.2f  |w_n|*sqrt(n)=%.3e  |I_n|=%.3e' % (t, a*math.sqrt(N)*t/(2*math.pi), abs(i)))
