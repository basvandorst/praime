"""Point 4 (final arithmetic, 50 digits) and the diagonal claim g = Psi(a,a)/2 on [1/1024, 7/8] from the referee's g."""
import os, numpy as np, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
JA = mp.mpf('1.6993342174476084346807187377464111453') + mp.mpf('1.28e-38')
cA = mp.mpf(3)/2*(2*JA*mp.mpf(9)/25)**(mp.mpf(1)/3)
D = mp.mpf(263470260835153127051917226381)/mp.mpf(397693215151001725278407884800)
JB = mp.mpf('1.3770253571246792139061286365308517464') + mp.mpf('2.75e-38')
cB = mp.mpf(3)/2*(2*JB*D*D)**(mp.mpf(1)/3)
print('A: c from certified J =', mp.nstr(cA, 20), ' <= 1.604335:', cA <= mp.mpf('1.604335'))
print('B: c from certified J =', mp.nstr(cB, 20), ' <= 1.597856:', cB <= mp.mpf('1.597856'))
# J = 4 pi (S0+S1+S2+S3) for A from the printed parts
S1 = mp.mpf('0.076188878925')
print('A: 4 pi (S0+S1+S2+S3) from printed parts =', mp.nstr(4*mp.pi*(mp.mpf('2.386514e-7')+mp.mpf('0.07618887892458584')+mp.mpf('0.058588735243')+mp.mpf('0.000450867517')), 15))
# referee float values
for name, J, DD in [('A', mp.mpf('1.6993340347'), mp.mpf(3)/5), ('B', mp.mpf('1.3770191667'), D)]:
    print(name, 'referee float c =', mp.nstr(mp.mpf(3)/2*(2*J*DD*DD)**(mp.mpf(1)/3), 12))
d = np.load(os.path.join(HERE, 'ref_g_ball.npz')); A, G = d['A'], d['G']
m = (A >= 1/1024) & (A <= 7/8)
a = A[m]; diag = a**6*(2-a)**4*(a*a+8*a+10)/320
print('ball: max over a in [1/1024,7/8] of (g - Psi(a,a)/2)/g =', np.max((G[m]-diag)/G[m]))
m2 = (A > 7/8) & (A < 0.95)
a = A[m2]; diag = a**6*(2-a)**4*(a*a+8*a+10)/320
dev = (G[m2]-diag)/G[m2]
print('ball: first a > 7/8 with g > Psi(a,a)/2 (rel 1e-12):', a[np.argmax(dev > 1e-12)])
