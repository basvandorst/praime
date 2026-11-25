"""CUE(N) simulation (Haar via QR of complex Ginibre with phase fix) of Y = log|Lambda'(theta_n)| against the exact
Palm-measure formulas of theory.py (E e^{tY} = M_{N-1}(1 + t/2)/M_{N-1}(1))."""
import numpy as np, sys
sys.path.insert(0, '.')
from theory import cue_cumulants, cue_moment
rng = np.random.default_rng(7)
for N, M in ((10, 30000), (20, 12000)):
    Y = []
    for _ in range(M):
        Zm = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2)
        Q, R = np.linalg.qr(Zm); d = np.diag(R); Q = Q * (d / np.abs(d))
        th = np.angle(np.linalg.eigvals(Q)); e = np.exp(1j * th)
        D = np.abs(e[:, None] - e[None, :]); np.fill_diagonal(D, 1.0); Y.append(np.log(D).sum(1))
    Y = np.concatenate(Y); c = Y - Y.mean(); m2 = (c ** 2).mean()
    sim = [Y.mean(), m2, (c ** 3).mean(), (c ** 4).mean() - 3 * m2 ** 2]; ex = cue_cumulants(N)
    se = [np.sqrt(m2 / len(Y)) * 3, np.sqrt(((c ** 2 - m2) ** 2).mean() / len(Y)) * 3]   # x3: eigenvalues of one matrix are correlated
    print('N=%d (%d matrices): mean %.4f (exact %.4f), var %.4f (%.4f), k3 %.4f (%.4f), k4 %.4f (%.4f)  [rough se %.4f, %.4f]' % (N, M, sim[0], ex[0], sim[1], ex[1], sim[2], ex[2], sim[3], ex[3], se[0], se[1]))
    for k in (-1, -0.5, 0.5, 1):
        print('   E|Lambda\'|^{2k}, k=%+.1f: sim %.5f exact %.5f' % (k, np.mean(np.exp(2 * k * Y)), cue_moment(N, k)))
