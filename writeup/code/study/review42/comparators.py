import numpy as np, sys, json
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from readings import *
import mypred
T = np.load(W + 'tables.npz')
P = mypred.Pred((T['lnh'], T['lnPhiD'], T['lnPhiO']))

def nodes(k, t0, t1, n=8):
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (t1 - t0) * x + 0.5 * (t1 + t0); w = w * 0.5 * (t1 - t0)
    Ls = np.log(k * t / (2 * np.pi)); wt = w * Ls
    return Ls, wt / wt.sum()

def law(ks, t0, t1, kind, dL=0.0):
    out = []
    for k in ks:
        bp = mypred.badprimes(k); Ls, wt = nodes(k, t0, t1)
        acc = 0
        for L, w in zip(Ls, wt):
            if kind == 'full':
                s = P.sigma2(lgrid, L + dL, bp, mode='full')
            elif kind == 'diag':
                s = P.sigma2(lgrid, L + dL, bp, mode='diag')
            elif kind == 'oneatatime':
                s = P.sigma2(lgrid, L, bp, mode='diag')
                for p in bp:
                    s = s + P.sigma2(lgrid, L, [p], mode='full') - P.sigma2(lgrid, L, [p], mode='diag')
            acc = acc + w * s
        out.append(acc)
    return np.array(out)

def gls_chi2(S, M, X):
    n, k = X.shape
    Y = (S - M)[:, pidx]
    b, *_ = np.linalg.lstsq(X, Y, rcond=None); e = Y - X @ b
    Sig = e.T @ e / (n - k); XtXi = np.linalg.inv(X.T @ X)
    V = np.kron(Sig, XtXi[1:, 1:]); r = b[1:].T.ravel()
    return float(r @ np.linalg.solve(V, r)), r, V

if __name__ == "__main__":
  for test in ("B", "A"):
      S, full, diag, X, primes, ks = load(test)
      t0, t1 = (100, 600) if test == 'B' else (100, 500)
      ks = [abs(k) for k in ks]
      M = {}
      M['full (BK)'] = full
      M['diag'] = diag
      M['sign-flipped off'] = 2 * diag - full
      M['one prime at a time'] = law(ks, t0, t1, 'oneatatime')
      for dL in (-0.3, +0.3, -0.6, +0.6):
          M['full at L%+.1f' % dL] = law(ks, t0, t1, 'full', dL)
      for a in (0.75, 0.85, 0.933, 1.15, 1.25):
          M['diag + %.3f*off' % a] = diag + a * (full - diag)
      print('=== test', test)
      o0 = readings(S, full, diag, X, primes)
      for nm, m in M.items():
          c2, r, V = gls_chi2(S, m, X)
          o = readings(S, m, diag, X, primes)   # registered statistics with m in place of full
          print('  %-22s GLS chi2 %7.1f on %d | OLS-se chi2 %7.1f | gamma(vs own off) %6.3f | rho pred %.3f obs %.3f | p2 l=6..8 pred %s' % (
              nm, c2, r.size, o['chi2f'], o['gamma'] if np.isfinite(o['gamma']) else np.nan, o['rho_pred'], o['rho'],
              ' '.join('%+.4f' % v for v in o['p2p'])))
      # (iv) diag + per-prime sinusoid of period 1 in l (a_p cos 2 pi l + b_p sin 2 pi l), GLS fit
      c2d, r, V = gls_chi2(S, diag, X)      # r = vec of (fit - pd), ordering l-major prime-minor
      npr = len(primes); L16 = prim
      A = np.zeros((16 * npr, 2 * npr))
      for j in range(16):
          for q in range(npr):
              A[j * npr + q, 2 * q] = np.cos(2 * np.pi * L16[j]); A[j * npr + q, 2 * q + 1] = np.sin(2 * np.pi * L16[j])
      Vi = np.linalg.inv(V)
      beta = np.linalg.solve(A.T @ Vi @ A, A.T @ Vi @ r)
      rr = r - A @ beta
      print('  diag + free sinusoid per prime (%d params): GLS chi2 %.1f on %d (diag alone %.1f)' % (2 * npr, rr @ Vi @ rr, r.size - 2 * npr, c2d))
      # same with a free per-prime smooth term too: a_p cos + b_p sin + c_p (constant) + d_p * l
      A2 = np.zeros((16 * npr, 4 * npr))
      for j in range(16):
          for q in range(npr):
              A2[j * npr + q, 4 * q:4 * q + 4] = [np.cos(2 * np.pi * L16[j]), np.sin(2 * np.pi * L16[j]), 1.0, L16[j]]
      beta2 = np.linalg.solve(A2.T @ Vi @ A2, A2.T @ Vi @ r); rr2 = r - A2 @ beta2
      print('  diag + sinusoid + linear per prime (%d params): GLS chi2 %.1f' % (4 * npr, rr2 @ Vi @ rr2))
      # BK off-shape with free amplitude per prime
      c2f, rf, Vf = gls_chi2(S, full, X)
      offv = (r - rf)   # = (pF - pD) vec
      A3 = np.zeros((16 * npr, npr))
      for j in range(16):
          for q in range(npr):
              A3[j * npr + q, q] = offv[j * npr + q]
      b3 = np.linalg.solve(A3.T @ Vi @ A3, A3.T @ Vi @ r); rr3 = r - A3 @ b3
      print('  diag + BK off-shape, free amplitude per prime: amplitudes %s, GLS chi2 %.1f (BK fixed: %.1f)' % (np.round(b3, 3), rr3 @ Vi @ rr3, c2f))
      # sinusoid amplitude estimates per prime vs BK
      print('  fitted sinusoid per prime (a cos + b sin):', [(p, round(beta[2 * q], 4), round(beta[2 * q + 1], 4)) for q, p in enumerate(primes)])
