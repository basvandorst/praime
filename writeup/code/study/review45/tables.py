# Tabulate regular parts of zeta(1+iy), zeta'/zeta(1+iy), (zeta'/zeta)'(1+iy) on y>0 (half-offset grid), mpmath.
import mpmath as mp, numpy as np, time
mp.mp.dps = 25
h = 0.002; Y = 6.6
J = int(Y/h)
ys = (np.arange(J) + 0.5)*h
gz = np.empty(J, complex); gl = np.empty(J, complex); gd = np.empty(J, complex)
t0 = time.time()
for j, y in enumerate(ys):
    s = mp.mpc(1, y)
    z0 = mp.zeta(s); z1 = mp.zeta(s, 0, 1); z2 = mp.zeta(s, 0, 2)
    x = mp.mpc(0, y)
    gz[j] = complex(z0 - 1/x)
    ld = z1/z0
    gl[j] = complex(ld + 1/x)
    gd[j] = complex(z2/z0 - ld**2 - 1/x**2)
print('time', time.time()-t0)
np.savez(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/ztab.npz', ys=ys, gz=gz, gl=gl, gd=gd)
print(gz[0], gl[0], gd[0])  # expect ~ euler gamma 0.5772, gamma 0.5772, ...
