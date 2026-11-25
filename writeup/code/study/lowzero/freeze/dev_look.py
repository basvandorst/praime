"""dev_look.py: first look at the development values: per character the means and per-interval spreads of the
observables, the discretisation of the maximum (grid vs parabola), and the lag-1 correlation between intervals."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obslib

DEVDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r149', 'dev_values')


def main():
    Ds = [1, 1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
    extra = sys.argv[1:]
    rows = []
    for D in Ds:
        v = np.load(os.path.join(DEVDIR, 'absL_%d.npy' % D))
        o = obslib.observables(v)
        q = abs(D)
        Lbar = np.mean(np.log(q * (100 + obslib.THETA * (np.arange(143) + 0.5)) / (2 * np.pi)))
        bad = [p for p in (2, 3, 5, 7) if D % p == 0 and D != 1]
        line = 'D=%6d L=%.2f bad=%-8s' % (D, Lbar, bad)
        for k in ['M', 'lnZ_0.5', 'lnZ_1', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'mean_ln']:
            x = o[k]
            r1 = np.corrcoef(x[:-1], x[1:])[0, 1]
            line += ' | %s %.3f sd %.3f r1 %+.2f' % (k, x.mean(), x.std(ddof=1), r1)
        dM = o['M'] - o['Mgrid']
        line += ' | dM mean %.4f max %.4f' % (dM.mean(), dM.max())
        print(line)


if __name__ == '__main__':
    main()
