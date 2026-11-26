# convergence of the predicted bin averages in the truncation parameters
import os, sys, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import predlib as PL
edges = np.round(np.linspace(0, 3, 61), 10)
print('#a values:', [(A, t, len(PL.enumerate_a(A, t))) for A, t in ((1e4, 1e-4), (1e5, 1e-5), (1e6, 1e-6))])
for sgn in (-1, 1):
    res = {}
    for (A, t, B0, Bc) in ((1e4, 1e-4, 200, 3000), (1e5, 1e-5, 300, 6000)):
        t0 = time.time(); m, tot = PL.predict_bins(edges, sgn, Amax=A, tau=t, B0=B0, Bcap=Bc, verbose=True)
        res[(A, t)] = m; print('sgn %+d Amax %g tau %g: %.1f s; bins 1..12:' % (sgn, A, t, time.time() - t0), np.round(m[:12], 3))
    a, b = list(res.values())
    print('   max |diff| over 60 bins: %.4f ; rms %.4f' % (np.max(np.abs(a - b)), np.sqrt(np.mean((a - b)**2))))
