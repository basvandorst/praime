# Point-mass problem (exploration after round 165): every prime carries ONE Satake angle (a genuine Euler product).
# t_pm(r) = max over angles (theta_p)_{p <= X} of min over forms of the smallest eigenvalue (per unit conductor exponent).
# Branch and bound: a node is a window of angles per prime; its upper bound is the convex problem with the measures
# restricted to the windows (cuts are shared between nodes); its lower bound is the value at one point of the box.
# Usage: pm.py <eps> <X> <r_lo> <r_hi> <step> [kmax] [maxnodes] [sign]     (eps = 2: holomorphic weight as parameter;
# sign: stop as soon as the sign of t_pm against -0.002 is decided)
import sys, json, os, time, heapq, numpy as np
from scipy.optimize import linprog
from maass import Problem

class PM:
    def __init__(self, r, eps, X, kmax=2, J=241, hol=False):
        self.P = P = Problem(r, eps, X, kmax=kmax, J=J, hol=hol)
        self.nP = len(P.primes); self.J = J; self.two_a = 2*P.B.a
        self.rows, self.rhs, self.tk = [], [], []
        for iT, A in enumerate(P.A):
            w, V = np.linalg.eigh(A)
            for i in range(min(V.shape[1], 12)): self.add(iT, V[:, i])
        self.Aeq = np.zeros((self.nP, self.nP*J + 1))
        for i in range(self.nP): self.Aeq[i, i*J:(i + 1)*J] = 1
    def add(self, iT, v):
        a_v, row = self.P.cut(iT, v); self.rhs.append(a_v); self.rows.append(row); self.tk.append(self.two_a*self.P.T[iT][2])
    def relax(self, win, tol=2e-5, rounds=200):
        P = self.P; J = self.J
        mask = np.zeros((self.nP, J))
        for i, (lo, hi) in enumerate(win): mask[i, lo:hi + 1] = 1
        cost = np.zeros(self.nP*J + 1); cost[-1] = -1
        bounds = [(0, float(m)) for m in mask.ravel()] + [(-50, 50)]
        for rd in range(rounds):
            res = linprog(cost, A_ub=np.hstack([np.array(self.rows), np.array(self.tk)[:, None]]), b_ub=np.array(self.rhs),
                          A_eq=self.Aeq, b_eq=np.ones(self.nP), bounds=bounds, method='highs')
            if res.status != 0: raise RuntimeError(res.message)
            w = res.x[:-1].reshape(self.nP, J); t = res.x[-1]; added = 0
            for iT, H in enumerate(P.forms(w)):
                ev, V = np.linalg.eigh(H); k = P.T[iT][2]
                for i in range(3):
                    if ev[i]/(self.two_a*k) < t - tol: self.add(iT, V[:, i]); added += 1
            if not added: break
        return t, w
    def point(self, th):
        """Exact value at angles th (one per prime); also adds the violated eigenvectors as cuts."""
        P = self.P; val = np.inf
        for iT, (A, x) in enumerate(zip(P.A, P.T)):
            cf = x[3]; coef = np.array([cf(m*th[P.primes.index(p)])*np.log(p)*p**(-m/2) for p, m in P.pp])
            ev, V = np.linalg.eigh(A - np.einsum('i,ijk->jk', coef, P.Spm))
            val = min(val, ev[0]/(self.two_a*x[2])); self.add(iT, V[:, 0])
        return val
    def solve(self, gap=2e-3, maxnodes=400, stop_at_sign=False):
        P = self.P; root = [(0, self.J - 1)]*self.nP
        t, w = self.relax(root); best = -np.inf; best_th = None; heap = [(-t, 0, root, w)]; cnt = 1; nodes = 1
        if stop_at_sign and t < -0.002: return -np.inf, t, np.array([float(w[i] @ P.th) for i in range(self.nP)]), 1      # the convex problem is infeasible already
        while heap:
            negt, _, win, w = heapq.heappop(heap); t = -negt
            if t <= best + gap: heap = []; ub = max(t, best); break
            th = np.array([float(w[i] @ P.th) for i in range(self.nP)])
            v = self.point(th)
            if v > best: best, best_th = v, th
            if stop_at_sign and (best >= -0.002 or t < -0.002): ub = t; break      # the sign of t_pm against -0.002 is decided
            if t <= best + gap: ub = t; break
            if nodes >= maxnodes: ub = t; break
            spread = [float(np.sqrt(max(0.0, w[i] @ (P.th - th[i])**2))) + 1e-3*(win[i][1] - win[i][0])/self.J for i in range(self.nP)]
            i = int(np.argmax([s if win[j][1] > win[j][0] else -1 for j, s in enumerate(spread)]))
            lo, hi = win[i]
            if hi == lo: ub = t; break
            mid = int(np.clip(np.searchsorted(P.th, th[i]), lo + 1, hi)) - 1
            for sub in ((lo, mid), (mid + 1, hi)):
                child = list(win); child[i] = sub
                tc, wc = self.relax(child); nodes += 1
                if tc > best + gap: cnt += 1; heapq.heappush(heap, (-tc, cnt, child, wc))
            ub = max([-h[0] for h in heap], default=best)
        else:
            ub = best
        return best, max(ub, best), best_th, nodes

if __name__ == '__main__':
    eps = int(sys.argv[1]); X = float(sys.argv[2]); lo, hi, step = map(float, sys.argv[3:6])
    kmax = int(sys.argv[6]) if len(sys.argv) > 6 else 2; maxnodes = int(sys.argv[7]) if len(sys.argv) > 7 else 400
    hol = eps == 2; sign = len(sys.argv) > 8 and sys.argv[8] == 'sign'
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass', 'pm_%s_X%g_k%d.jsonl' % ('hol' if hol else 'eps%d' % eps, X, kmax))
    done = set()
    if os.path.exists(out):
        for line in open(out): done.add(round(json.loads(line)['r'], 4))
    for r in [round(lo + i*step, 4) for i in range(int(round((hi - lo)/step)) + 1)]:
        if r in done: continue
        t0 = time.time(); S = PM(r, 0 if hol else eps, X, kmax=kmax, hol=hol); best, ub, th, nodes = S.solve(maxnodes=maxnodes, stop_at_sign=sign)
        rec = dict(r=r, eps=eps, X=X, kmax=kmax, mode='sign' if sign else 'full', t_point=best if np.isfinite(best) else None, t_upper=ub, nodes=nodes, secs=round(time.time() - t0, 1), primes=S.P.primes,
                   ap=[float(2*np.cos(x)) for x in th], cuts=len(S.rows))
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
        print('r %.4f: point-mass t in [%+.5f, %+.5f], %d nodes, %.0fs, a_p %s' % (r, best, ub, nodes, rec['secs'], ' '.join('%d:%+.3f' % (p, a) for p, a in zip(S.P.primes, rec['ap']))), flush=True)
