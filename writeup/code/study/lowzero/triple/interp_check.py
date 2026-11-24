"""check of the cached/interpolated path of predict3.kf3_parts against the direct tplib computation at off-grid L"""
import sys, numpy as np
src = open('predict3.py').read().split("out = {'l': LS")[0].replace("MODE = sys.argv[1]\nif MODE == 'tables': sys.exit(0)", "MODE='none'")
exec(src)
PS = tplib.PrimeSums(10 ** 4); worst = 0
for L, S in ((10.039, (2, 3, 89)), (10.762, (3, 373)), (10.112, (3, 383)), (9.437, ()), (8.71, (5,)), (11.546, (2, 7)), (7.93, (3,))):
    Qp, Pp = kf3_parts(L, S); d = tplib.T3_parts(0 * q.T, q.T, q.S, L, S, PS=PS); Pd = q.kf3(d['P']); Qd = q.kf3(d['Q'])
    e = max(np.abs(Pp - Pd).max(), np.abs(Qp - Qd).max()); worst = max(worst, e)
    print('L=%.3f S=%-12s max |cached - direct| over l: P %.1e  Q %.1e' % (L, S, np.abs(Pp - Pd).max(), np.abs(Qp - Qd).max()))
print('largest %.1e' % worst)
