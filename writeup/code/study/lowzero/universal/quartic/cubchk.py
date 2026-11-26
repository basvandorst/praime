# cubic check of the Mellin evaluation against the round-161 table (Poisson-summed closed forms), bin by bin
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P
tab = np.loadtxt('/home/bas/praime/writeup/code/study/lowzero/universal/cubic/m_table.txt')
cubs = {T: P.cubic_model(T=T) for T in (1000.0, 2000.0)}
for sg, col in ((-1, 2), (1, 3)):
    # bins of width 0.03 with edges at 0.01 + 0.03 j (away from the spikes k^2/4, k^2/8 where possible)
    e = np.round(np.arange(0.01, 2.99, 0.03), 10); ref = np.array([tab[int(round(a * 100)):int(round(b * 100)), col].mean() for a, b in zip(e[:-1], e[1:])])
    for T, cm in cubs.items():
        mine = P.bin_means(cm, P.CUBIC_MU[sg], e); d = mine - ref; k = np.argsort(-np.abs(d))[:4]
        print('D%s0 T=%5.0f: bins of 0.03: max |diff| %.4f rms %.4f; worst at y = %s' % ('<' if sg < 0 else '>', T, np.abs(d).max(), np.sqrt(np.mean(d**2)), [round(float(e[i]), 2) for i in k]))
