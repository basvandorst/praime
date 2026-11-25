import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, numpy as np
from smlib import *
N, k = 10007, 4
for T in [float(a) for a in sys.argv[1:]]:
    win = Window(N, k, T, h=0.05)
    Hh = lambda xi: np.sum(win.H*np.exp(-1j*xi*win.tt))*win.h
    ms = np.arange(1, 400)
    v = np.array([Hh(-math.log(m)) for m in ms])
    p11 = np.sum(np.abs(v)**2/ms)/(4*math.pi**2); p2 = np.sum(v**2/ms)/(4*math.pi**2)
    print('T=%g: integer model theorem: M11 = %.5f, M2 = %.5f%+.5fi (L^2 = %.5f; m=2 term |H^(ln2)|/H^(0) = %.3f)' % (T, p11, p2.real, p2.imag, win.L**2, abs(v[1])/abs(v[0])))
