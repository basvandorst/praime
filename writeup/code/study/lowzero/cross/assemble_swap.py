"""Assemble the frozen swapped-term predictions (primes <= 5e4, u-step 0.06) into swap_pred_NN.npz, swap_pred_NO.npz;
print their sha256, the number of pairs, the sum of squares of the predictor and its size distribution."""
import os, sys, glob, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for des in ('NN', 'NO'):
    fs = sorted(glob.glob(os.path.join(HERE, 'swap_%s_P50000_*.npz' % des)))
    D1 = np.concatenate([np.load(f)['D1'] for f in fs]); D2 = np.concatenate([np.load(f)['D2'] for f in fs])
    F = np.concatenate([np.load(f)['F'] for f in fs]); C = np.concatenate([np.load(f)['C'] for f in fs])
    out = os.path.join(HERE, 'swap_pred_%s.npz' % des); np.savez(out, D1=D1, D2=D2, F=F, C=C)
    sha = hashlib.sha256(open(out, 'rb').read()).hexdigest()
    a = np.abs(C)
    print('%s: %d pairs (F > 0.1), sum of squares %.4f, |C| median %.4f, 90%% %.4f, max %.4f; sha256 %s' % (des, len(D1), float((C ** 2).sum()), np.median(a), np.percentile(a, 90), a.max(), sha))
    print('   per w rms: ' + ' '.join('%.4f' % v for v in np.sqrt((C ** 2).mean(0))))
