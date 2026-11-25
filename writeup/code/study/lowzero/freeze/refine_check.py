"""refine_check.py: discretisation of the interval maximum.  For development characters, at every interval's grid
maximum t_km (step 0.05), evaluate |L(1/2+it)| with PARI on 41 points in [t_km - 0.05, t_km + 0.05] (step 0.0025),
take the maximum there refined by a parabola (error < 1e-6), and compare with (a) the plain grid maximum and (b) the
parabola through the three grid points used by obslib/freeze_test.  usage: refine_check.py D [D ...]"""
import sys, os
import numpy as np
import cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obslib

DEV8 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r149', 'dev_values')
DEVX = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_extra_values')


def main():
    pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6); pari.set_real_precision(19)
    for D in [int(x) for x in sys.argv[1:]]:
        f = os.path.join(DEV8, 'absL_%d.npy' % D)
        if not os.path.exists(f):
            f = os.path.join(DEVX, 'absL_%d.npy' % D)
        v = np.load(f)
        o = obslib.observables(v)
        pari('Li = lfuninit(lfuncreate(%d), [1/2, 0, 1001])' % D)
        g = obslib.grid(); lnL = np.log(v)
        d_grid, d_par = [], []
        for j in range(obslib.NINT):
            a = obslib.T0 + j * obslib.THETA
            k0 = int(np.ceil((a - obslib.T0) / obslib.H - 1e-9))
            k1 = min(int(np.floor((a + obslib.THETA - obslib.T0) / obslib.H + 1e-9)), obslib.NPT - 1)
            km = k0 + int(np.argmax(lnL[k0:k1 + 1]))
            if km in (k0, k1):
                continue                      # maximum at an interval end: skip (parabola not used there)
            tt = g[km] + 0.0025 * np.arange(-20, 21)
            vals = np.array([float(pari('log(abs(lfun(Li, 1/2 + I*%.10f)))' % x)) for x in tt])
            i = int(np.argmax(vals))
            if 0 < i < 40:
                y0, y1, y2 = vals[i - 1:i + 2]
                den = y0 - 2 * y1 + y2
                true = y1 - (y0 - y2) ** 2 / (8 * den) if den < 0 else y1
            else:
                true = vals[i]
            d_grid.append(true - lnL[km]); d_par.append(true - o['M'][j])
        d_grid = np.array(d_grid); d_par = np.array(d_par)
        print('D=%6d  intervals %d | true - grid max: mean %.5f max %.5f | true - parabola: mean %+.6f '
              'sd %.6f max|.| %.6f' % (D, len(d_grid), d_grid.mean(), d_grid.max(), d_par.mean(), d_par.std(),
                                       np.abs(d_par).max()), flush=True)


if __name__ == '__main__':
    main()
