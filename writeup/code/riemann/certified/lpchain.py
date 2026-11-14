"""LP-certified chain: primes proved from zeta zeros, decisions certified by Farkas vectors.

Same samples as chain.py: y(x) = sum_{u in U} b_u ln(u) sqrt(x/u) K(x ln(u/x)) + e(x), |e(x)| <= delta(x),
b in {0,1}^k. A value v of unknown j is certified when b_j = 1 - v is refuted:
branch and bound over the other unknowns; every leaf is closed by
  (i)  propagation: all coefficients are >= 0, so a row forces b_l = 0 or 1 or is violated outright, or
  (ii) a Farkas vector w over the samples: min over the box of (A^T w).b  >  w.y + |w|.delta,
       which no b in the box can satisfy (checked in floating point with a margin of 1e-9 relative).
If a leaf holds an integer pattern that satisfies every sample, the value is ambiguous: the run stops.
Exploration mode (truth_below = S) takes the true Lambda below S; the chain from 3 uses no prime as input.
"""
import math, sys, time, json, os
import numpy as np
from scipy.optimize import linprog
from kbm import K, K0
from cert import sample

HERE = os.path.dirname(os.path.abspath(__file__))
ZFILE = os.environ.get('ZFILE', os.path.join(HERE, '../../../data/riemann/zeros_hp2419.txt'))
if ZFILE.endswith('.txt'):      # decimal strings, read into 64-bit-mantissa floats (error < 1.1e-16 here)
    G = np.array([l.strip() for l in open(ZFILE) if l.strip()], dtype=np.longdouble)
else:
    G = np.load(ZFILE) if ZFILE.endswith('.npy') else np.fromfile(ZFILE)
EPS = float(os.environ.get('EPS', '2e-16' if G.dtype == np.longdouble else '1e-12'))
REL = float(os.environ.get("REL", "1e-12"))   # margin for floating-point comparisons that close a leaf (inputs carry < 1e-13 relative error; checked in mpmath)

def iroot(n, e):
    r = int(round(n ** (1.0 / e)))
    for q in (r - 1, r, r + 1):
        if q > 1 and q ** e == n:
            return q
    return 0

def propagate(A, lo, hi, st, rec=None):
    """Fix variables forced by single rows (A >= 0). Returns False if a row is violated.
    rec (optional list) receives ('force', var, val, row) and ('viol', row, side) for the checker."""
    while True:
        one = st == 1; free = np.where(st < 0)[0]
        base = A[:, one].sum(1)
        Af = A[:, free]; tot = Af.sum(1)
        tol = REL * (1 + np.abs(hi) + np.abs(lo) + base + tot)
        bad_hi = base > hi + tol; bad_lo = base + tot < lo - tol
        if bad_hi.any() or bad_lo.any():
            if rec is not None:
                rec.append(('viol', int(np.argmax(bad_hi)), 'hi') if bad_hi.any() else ('viol', int(np.argmax(bad_lo)), 'lo'))
            return False
        if len(free) == 0:
            return True
        O = base[:, None] + Af > (hi + tol)[:, None]
        Un = base[:, None] + tot[:, None] - Af < (lo - tol)[:, None]
        over = O.any(0); under = Un.any(0)
        if np.any(over & under):
            j = int(np.argmax(over & under))
            if rec is not None:
                rec.append(('force', int(free[j]), 0, int(np.argmax(O[:, j]))))
                rec.append(('viol', int(np.argmax(Un[:, j])), 'lo1'))
            return False
        if not over.any() and not under.any():
            return True
        if rec is not None:
            for j in np.where(over)[0]: rec.append(('force', int(free[j]), 0, int(np.argmax(O[:, j]))))
            for j in np.where(under)[0]: rec.append(('force', int(free[j]), 1, int(np.argmax(Un[:, j]))))
        st[free[over]] = 0; st[free[under]] = 1

def farkas(Af, lo_r, hi_r):
    """LP over w: maximise min_box (Af^T w).b - w.y - |w|.d with ||w||_1 = 1.
    Returns (True, w) if the verified value is positive (infeasible box), else (False, b*) with b* the
    LP point of least violation (the dual solution)."""
    keep = Af.max(1) > 0
    Af, lo_r, hi_r = Af[keep], lo_r[keep], hi_r[keep]
    n, kf = Af.shape
    farkas.keep = keep
    yv = (lo_r + hi_r) / 2; dv = (hi_r - lo_r) / 2
    c = np.concatenate([yv + dv, -yv + dv, -np.ones(kf)])
    A_ub = np.hstack([-Af.T, Af.T, np.eye(kf)])
    A_eq = np.concatenate([np.ones(2 * n), np.zeros(kf)])[None]
    bounds = [(0, None)] * (2 * n) + [(None, 0)] * kf
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(kf), A_eq=A_eq, b_eq=[1.0], bounds=bounds, method='highs')
    if res.status != 0:
        return False, np.full(kf, 0.5)
    w = res.x[:n] - res.x[n:2 * n]
    g = Af.T @ w
    minbox = np.minimum(g, 0).sum()
    rhs = w @ yv + np.abs(w) @ dv
    scale = np.abs(w) @ (Af.sum(1) + np.abs(yv) + dv)
    if minbox - rhs > REL * scale:
        return True, w
    b = -np.asarray(res.ineqlin.marginals)
    return False, np.clip(b, 0, 1)

def refute(A, lo, hi, st0, limit, log=None):
    """Branch and bound: is there b in {0,1}^k, b = st0 where st0 >= 0, with lo <= A b <= hi?
    Returns ('refuted', nodes, lps), ('feasible', pattern), or ('limit', nodes, lps).
    log (optional list) receives one record per node: the state before propagation, the forcings,
    and how the node was closed or split; a separate checker replays it."""
    stack = [st0.copy()]; nodes = 0; lps = 0
    while stack:
        st = stack.pop(); nodes += 1
        if nodes > limit:
            return ('limit', nodes, lps)
        rec = None if log is None else []
        node = None if log is None else dict(state=''.join('x01'[v + 1] for v in st), events=rec)
        if node is not None: log.append(node)
        if not propagate(A, lo, hi, st, rec):
            continue
        free = np.where(st < 0)[0]
        if len(free) == 0:
            return ('feasible', st.copy())
        base = A[:, st == 1].sum(1)
        lps += 1
        ok, v = farkas(A[:, free], lo - base, hi - base)
        if ok:
            if node is not None:
                w = np.zeros(A.shape[0]); w[np.where(farkas.keep)[0]] = v; node['farkas'] = w.tolist()
            continue
        j = 0                      # branch in index order: the samples are nearly triangular
        r = int(round(v[j]))
        if node is not None: node['branch'] = int(free[j])
        for val in (1 - r, r):     # the LP-rounded child is explored first
            s2 = st.copy(); s2[free[j]] = val; stack.append(s2)
    return ('refuted', nodes, lps)

class LPChain:
    def __init__(self, X, wheel=(2,), truth_below=0):
        self.X = X; N = X + 2000
        self.wheel = wheel
        self.lam = np.full(N, np.nan)
        self.lam[0] = self.lam[1] = 0.0
        W = 1
        for p in wheel: W *= p
        self.W = W
        for n in range(2, N):
            for p in wheel:
                if n % p == 0:
                    q = n
                    while q % p == 0: q //= p
                    self.lam[n] = math.log(p) if q == 1 else 0.0
                    break
        self.powers = {}
        for e in range(2, 40):
            r = 2
            while r ** e < N:
                if math.gcd(r, W) == 1 and not any(iroot(r, f) for f in range(2, 20)):
                    self.powers.setdefault(r, []).append(r ** e)
                r += 1
        self.is_pp = np.zeros(N, dtype=bool)
        for r, lst in self.powers.items():
            for v in lst: self.is_pp[v] = True
        if truth_below:
            from sympy import isprime
            for n in range(2, truth_below):
                if np.isnan(self.lam[n]) and not self.is_pp[n]:
                    self.commit(n, isprime(n))
        self.cache = {}; self.maxT = 0.0; self.log = []; self.stats = dict(nodes=0, lps=0, maxnodes=0)
        self.record = (); self.proofs = {}; self.cur = None; self.dlog = []

    def commit(self, n, prime):
        self.lam[n] = math.log(n) if prime else 0.0
        for v in self.powers.get(n, []):
            self.lam[v] = math.log(n) if prime else 0.0

    def zs(self, x, a, beta, m, c):
        key = (x, a, beta, m, c)
        v = self.cache.get(key)
        if v is None:
            T = c * beta * x / a
            self.maxT = max(self.maxT, T)
            z, d, _ = sample(x, G, T, a, beta, m, EPS)
            v = self.cache[key] = (z, d)
        return v

    def unknowns(self, c, k):
        U = []; n = c
        while len(U) < k + 1:
            if np.isnan(self.lam[n]) and not self.is_pp[n]:
                U.append(n)
            n += 1
        return U[:-1], U[-1]

    def window(self, c, a, beta, m, cc, k, h):
        U, w = self.unknowns(c, k)
        Ua = np.array(U, dtype=float)
        xs = []
        x = max(2.0, math.floor(c - a - 1))
        while x < w:
            if x * math.exp(-a / x) > 1.0 and abs(x * math.log(c / x)) < a and x * math.log(w / x) >= a:
                xs.append(float(x))
            x += h
        k0 = K0(beta, m)
        A = np.empty((len(xs), len(U))); y = np.empty(len(xs)); d = np.empty(len(xs))
        for i, x in enumerate(xs):
            z, dd = self.zs(x, a, beta, m, cc)
            lo, hi = math.floor(x * math.exp(-a / x)), math.ceil(x * math.exp(a / x))
            mm = np.arange(max(2, lo), hi + 1)
            v = x * np.log(mm / x)
            inside = np.abs(v) < a
            lv = self.lam[mm]
            if np.any(np.isnan(lv) & inside & ~np.isin(mm, Ua)):
                raise RuntimeError(f'unknown number inside the support of sample {x}')
            ok = inside & ~np.isnan(lv) & (lv > 0)
            known = float(np.sum(lv[ok] * np.sqrt(x / mm[ok]) * K(v[ok], a, beta, m)))
            A[i] = np.log(Ua) * np.sqrt(x / Ua) * K(x * np.log(Ua / x), a, beta, m) / k0
            y[i] = (z - known) / k0; d[i] = dd / k0
        self.last_xs = xs
        return U, A, y - d, y + d

    def decide(self, U, A, lo, hi, limit, maxbits):
        """Certify the leading unknowns one by one: the LP point suggests a value v for unknown j,
        and b_j = 1 - v is refuted (or, failing that, b_j = v). Returns the certified values."""
        k = len(U)
        logs = [] if self.record and self.record[0] <= U[0] < self.record[1] else None
        got = []; st = np.full(k, -1)
        for j in range(min(k, maxbits)):
            s0 = st.copy()
            if not propagate(A, lo, hi, s0):
                return got, 'data inconsistent'
            if s0[j] >= 0:
                guess = int(s0[j])
            else:
                free = np.where(s0 < 0)[0]
                base = A[:, s0 == 1].sum(1)
                ok, v = farkas(A[:, free], lo - base, hi - base)
                if ok:
                    return got, 'data inconsistent'
                guess = int(round(v[int(np.where(free == j)[0][0])]))
            done = False; outcomes = []; attempts = []; fits = {}
            for val in (guess, 1 - guess):
                s2 = st.copy(); s2[j] = 1 - val
                lg = None if logs is None else []
                rr = refute(A, lo, hi, s2, limit, lg)
                attempts.append((1 - val, rr[0], int(rr[1]) if rr[0] != 'feasible' else -1))
                if rr[0] == 'feasible':
                    fits[1 - val] = [int(v) for v in rr[1]]      # a 0/1 pattern with b_j = 1 - val that fits every sample
                if rr[0] != 'feasible':
                    self.stats['nodes'] += rr[1]; self.stats['lps'] += rr[2]
                    self.stats['maxnodes'] = max(self.stats['maxnodes'], rr[1])
                if rr[0] == 'refuted':
                    if logs is not None:
                        logs.append(dict(bit=j, value=val, root=''.join('x01'[v + 1] for v in s2), nodes=lg))
                    got.append(val); st[j] = val; done = True
                    break
                outcomes.append(rr[0])   # 'feasible': a pattern with b_j = 1 - val fits; try the other value
            self.dlog.append(dict(u=int(U[j]), j=j, c=int(U[0]), k=k, guess=int(guess), value=(got[-1] if done else None), attempts=attempts, a=self.cur[0] if self.cur else None))
            if not done:
                if logs is not None and got:
                    self.proofs[U[0]] = dict(U=U, xs=self.last_xs, params=self.cur, proofs=logs)
                if len(fits) == 2:   # the information limit: two patterns that differ in b_j both fit (checked in info_law_check.py)
                    self.ambiguous = dict(U=[int(u) for u in U], xs=self.last_xs, params=self.cur, bit=j, fits={str(v): p for v, p in fits.items()})
                return got, ('ambiguous' if outcomes == ['feasible', 'feasible'] else 'search limit ' + str(outcomes))
        if logs is not None:
            self.proofs[U[0]] = dict(U=U, xs=self.last_xs, params=self.cur, proofs=logs)
        return got, None

    def run(self, stages, start=3, limit=400, maxbits=8, h=1.0, verbose=False, record=()):
        c = start; t0 = time.time(); decided = 0; self.record = tuple(record)
        for (end, a, beta, m, cc, k) in stages:
            while c < min(end, self.X):
                if not np.isnan(self.lam[c]) or self.is_pp[c]:
                    if self.is_pp[c] and np.isnan(self.lam[c]):
                        raise RuntimeError(f'perfect power {c} reached before its root was decided')
                    c += 1; continue
                self.cur = [a, beta, m, cc, k, h]
                U, A, lo, hi = self.window(c, a, beta, m, cc, k, h)
                got, why = self.decide(U, A, lo, hi, limit, maxbits)
                if not got:
                    return dict(stop=c, reason=why, decided=decided, secs=round(time.time() - t0, 1), **self.stats)
                for u, v in zip(U, got):
                    self.commit(u, v == 1); decided += 1
                c = U[len(got) - 1] + 1
                if verbose and decided % 200 < len(got):
                    print(c, decided, round(time.time() - t0), self.stats, flush=True)
            self.log.append(dict(end=end, a=a, beta=beta, m=m, c=cc, k=k, height_over_pi_x=cc * beta / a / math.pi))
        return dict(stop=None, decided=decided, maxT=self.maxT, secs=round(time.time() - t0, 1), **self.stats)

if __name__ == '__main__':
    X = int(sys.argv[1]); stages = json.loads(sys.argv[2])
    opts = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    wheel = tuple(opts.pop('wheel', [2])); tb = opts.pop('truth_below', 0); start = opts.pop('start', 3)
    out = opts.pop('out', None)
    C = LPChain(X, wheel, tb)
    r = C.run(stages, start=start, **opts)
    print(r)
    from sympy import isprime
    lo = max(start, 3)
    primes = [n for n in range(lo, X) if C.lam[n] == math.log(n)]
    hi_dec = r['stop'] or X
    wrong = [n for n in range(lo, hi_dec) if not np.isnan(C.lam[n]) and (C.lam[n] == math.log(n)) != isprime(n)]
    print(f'{len(primes)} primes certified in [{lo}, {hi_dec}); wrong: {wrong[:10]}; zeros up to {C.maxT:.1f}')
    if out:
        json.dump(dict(X=X, stages=stages, wheel=wheel, start=start, result=r, log=C.log, primes=primes, opts=opts), open(out, 'w'))
        json.dump(C.dlog, open(out.replace('.json', '_decisions.json'), 'w'))
        if getattr(C, 'ambiguous', None):
            json.dump(C.ambiguous, open(out.replace('.json', '_ambiguous.json'), 'w'))
        if C.proofs:
            json.dump({str(k): v for k, v in C.proofs.items()}, open(out.replace('.json', '_proofs.json'), 'w'))
