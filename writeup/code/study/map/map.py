#!/usr/bin/env python3
"""map.py -- direction map of the within part of transition tables of primes mod q, with a linear-response comparator.

Self-contained: numpy and the standard library only.  Nothing is imported from the repo.

THE STATISTIC (exactly the forty-fourth review's dirs.py, generalised to any q and lag k)
  Window [x0, x0 + 1e8), blocks of length H (1e4, 1e5, 1e6) by the position of p_n.  c_n = class of p_n mod q.
  X_ab,b = #{n : p_n in block b, c_n = a, c_(n+k) = b}  (table of the block), T = sum over blocks (window table),
  P^ = T row-normalised (per window), E = N_a,b P^_ab, y_b = (X - E)/sqrt(E) in R^(C^2).
  Basis: for every pair of non-principal characters (chi_i, chi_j) of (Z/q)*, the function chi_i(a) conj chi_j(b),
  weighted by sqrt(pi_a P_ab) and projected on the within space (row and column parts removed, Pi_W of
  within/derivation.md 2.1); real and imaginary parts; Loewdin-orthonormalised.  P and pi are pooled over all windows
  of the set.  (C - 1)^2 real directions d.
  W_d(H) per window = sum_b (y_b . B_d)^2 / (NB - 1)   (= mean over blocks times NB/(NB - 1)).
  Measured slope per window and direction: s_d = (W_d(1e6) - W_d(1e4)) / ln 100.
  Per window we store M = sum_b y_b y_b^T (C^2 x C^2), so any basis can be applied afterwards.

THE COMPARATOR (linear response, the review's tailemp.py, per direction)
  Inserting a prime at a site x (coprime to 510510) changes the table by r_x = sum_added F(pair) - sum_removed F(pair),
  F(a, b) = within coordinates of the single-transition increment ((e_ab - P_a)/sqrt(pi_a P_ab)).  k = 1: added
  (L1, x), (x, R1), removed (L1, R1); general k: added (L_i, R_(k-i)) i < k, (L_k, x), (x, R_k); removed
  (L_i, R_(k+1-i)).  r_u = mean of r_x over sites x = u mod Q, Q = lcm(prod_(p<=y) p, q), y = 17 (11, 13 as checks).
  A_d(h) = (1/phi(Q)) sum_u r_(u,d) r_(u+h,d) from the class means pooled over the whole response stretch, except at
  h = 0 (and h = Q), where the sampling noise is removed by using cross products of the means of 4 disjoint contiguous
  quarters of the stretch instead (the review used var/(n-1); the two agree to 0.04 %).  Errors of the comparator:
  jackknife over the 4 quarters.  A check estimator ('cross', cross products of quarters at every h) is also stored.
  W_tail,d(H) = (2/lam) (Q/phi(Q)) (lam_resp/lam) sum_(h<H) (1 - h/H) (1 - e^(-h/lam)) A_d(h mod Q) (S_(>y)(h) - 1)
  (S_(>y) = Hardy-Littlewood pair factor of the primes > y; lam_resp/lam scales V ~ 1/lam to the measured height).
  LR slope = (W_tail,d(1e6) - W_tail,d(1e4)) / ln 100;  asymptotic variant -V_d/lam, V_d = A_d(0) (dirs.py).
  For k = 1 and C <= 6 the review's own noise-corrected estimate (pooled counts, A(0) -= mean var/(n-1)) is also given.
  The comparator at y is skipped when q has a prime factor above y (q = 13 at y = 11).
  EXCESS per direction d = mean over windows of s_d - LR slope_d; per type: mean over the directions of the type.
  Errors: jackknife over windows (measurement) and over the 4 quarters of the response stretch (comparator), combined
  in quadrature.

CONTROLS  measure --control perm: successors permuted among the pairs with the same first class in the whole window
  (expected W = 1 and slope 0 in every direction: checks the normalisation of every direction type).
  measure --control bshuf: the 1e4-blocks of each window in random order before they are summed to 1e5 and 1e6
  (E[slope] = 0 exactly for any data: an empirical null for the error bars).
  chain: a first-order chain with the primes' pooled matrix (no long range: the slope is the O(1/N) artefact).
  cells q:50 (p_(n+50) instead of p_(n+1)).

DIRECTION TYPES  eq: i = j;  conj: i = conj j (i != j);  rr: both real, different;  rc: one real, one complex;
  cc: both complex, neither equal nor conjugate.  Subtype: (ord chi_i, ord chi_j, ord chi_i conj chi_j).

USAGE (each run under: (ulimit -v 1500000; timeout 600 nice python -I map.py ...))
  map.py measure  OUTDIR x0 [x0 ...] [--cells SPEC] [--control none|perm|bshuf] [--seed S]
  map.py respond  --out LR.json --stretch x0 [x0 ...] --set NAME=DIR[,DIR] [--set ...] --qs 5 --ks 1,2
                  [--raw-out RAW.npz] [--raw-in RAW.npz ...] [--parts 4]   (several stretches are pooled)
  map.py chain    OUTDIR --set NAME=DIR[,DIR] --nwin 24 --seed S [--cells SPEC] [--first I]
  map.py report   --set NAME=DIR[,DIR] [--lr LR.json ...] [--json OUT.json] [--groups]
  map.py combine  REPORT1.json REPORT2.json [--json OUT]     (inverse-variance pooling of two heights)
  SPEC: comma list of q:k (default all q in 3,5,7,8,9,11,12,13 with k = 1, 2).
  Time and memory (measured): measure 2-5 s per window (all 18 cells); respond per q and 4 stretches <= 160 s and
  <= 1.1 GB resident (q = 13); report < 60 s.
"""
import sys, os, math, json, argparse, itertools, cmath, time, glob
import numpy as np

SPAN = 10 ** 8
PAD = 20000                      # enough for 50 successors at these heights
HS = (10 ** 4, 10 ** 5, 10 ** 6)
QS_ALL = (3, 5, 7, 8, 9, 11, 12, 13)
GENS = {3: [(2, 2)], 5: [(2, 4)], 7: [(3, 6)], 8: [(3, 2), (5, 2)], 9: [(2, 6)], 10: [(3, 4)],
        11: [(2, 10)], 12: [(5, 2), (7, 2)], 13: [(2, 12)]}
YQ = {11: 2310, 13: 30030, 17: 510510}
LN100 = math.log(100.0)


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------ primes
def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)


def sieve(lo, hi, base):
    """primes in [lo, hi), odd-only segmented sieve (base must hold the primes up to sqrt(hi))"""
    out = []
    SEG = 2 * 10 ** 7
    for a in range(lo, hi, SEG):
        b = min(hi, a + SEG)
        a0 = a | 1
        n = (b - a0 + 1) // 2
        if n <= 0: continue
        seg = np.ones(n, bool)
        for p in base[1:].tolist():
            if p * p >= b: break
            st = max(p * p, ((a0 + p - 1) // p) * p)
            if st % 2 == 0: st += p
            seg[(st - a0) // 2::p] = False
        out.append(a0 + 2 * np.flatnonzero(seg))
    P = np.concatenate(out)
    return P[P >= max(lo, 3)]


def primes_around(x0):
    base = small_primes(int(math.isqrt(x0 + SPAN + PAD)) + 2)
    P = sieve(x0 - PAD, x0 + SPAN + PAD, base)
    assert P[0] < x0 and P[-1] >= x0 + SPAN
    return P


# ------------------------------------------------------------------ the group and its characters
class Group:
    def __init__(self, q):
        self.q = q
        self.red = [a for a in range(q) if math.gcd(a, q) == 1]; C = self.C = len(self.red)
        gens = GENS[q]; ords = [o for _, o in gens]
        lg = {}
        for ex in itertools.product(*[range(o) for o in ords]):
            x = 1
            for (g, _), e in zip(gens, ex): x = x * pow(g, e, q) % q
            lg[x] = ex
        assert sorted(lg) == self.red, (q, sorted(lg))
        self.keys = keys = list(itertools.product(*[range(o) for o in ords]))
        self.tab = np.array([[cmath.exp(2j * math.pi * sum(kk * e / o for kk, e, o in zip(key, lg[a], ords)))
                              for a in self.red] for key in keys])
        self.conj = [keys.index(tuple((-kk) % o for kk, o in zip(key, ords))) for key in keys]
        self.order = [math.lcm(*[o // math.gcd(kk, o) for kk, o in zip(key, ords)]) for key in keys]
        self.real = [self.conj[i] == i for i in range(C)]
        self.ords = ords
        self.lut = -np.ones(q, np.int64); self.lut[self.red] = np.arange(C)
        # conductor of each character
        self.cond = []
        for i in range(C):
            for d in sorted(d for d in range(1, q + 1) if q % d == 0):
                if all(abs(self.tab[i][self.red.index(a)] - 1) < 1e-9 for a in self.red if a % d == 1 % d):
                    self.cond.append(d); break
        # orthogonality check
        G = self.tab @ self.tab.conj().T
        assert np.allclose(G, C * np.eye(C))
        self.groups = self._groups()

    def ratio(self, i, j):
        ki, kj = self.keys[i], self.keys[j]
        return self.keys.index(tuple((a - b) % o for a, b, o in zip(ki, kj, self.ords)))

    def name(self, i):
        return '%d' % self.keys[i][0] if len(self.keys[i]) == 1 else ''.join(str(k) for k in self.keys[i])

    def _groups(self):
        C = self.C; done = set(); out = []
        cyc = len(self.ords) == 1
        for i in range(1, C):
            for j in range(1, C):
                if (i, j) in done: continue
                ic, jc = self.conj[i], self.conj[j]; done.add((i, j)); done.add((ic, jc))
                n = 1 if (i, j) == (ic, jc) else 2
                if i == j: typ = 'eq'
                elif i == self.conj[j]: typ = 'conj'
                elif self.real[i] and self.real[j]: typ = 'rr'
                elif self.real[i] or self.real[j]: typ = 'rc'
                else: typ = 'cc'
                sub = (self.order[i], self.order[j], self.order[self.ratio(i, j)])
                if cyc:   # the review's label
                    rv = 'ratio' if i == j else ('product' if (i + j) % C == 0 else 'mixed')
                    label = '(%d,%d)%s' % (i, j, rv)
                else:
                    label = '(%s,%s)' % (self.name(i), self.name(j))
                out.append(dict(i=i, j=j, n=n, type=typ, sub='%d.%d.%d' % sub, label=label,
                                cond='%d,%d' % (self.cond[i], self.cond[j])))
        assert sum(g['n'] for g in out) == (C - 1) ** 2
        return out


def maxpf(q):
    return max(p for p in range(2, q + 1) if q % p == 0 and all(p % r for r in range(2, p)))


_GROUPS = {}
def group(q):
    if q not in _GROUPS: _GROUPS[q] = Group(q)
    return _GROUPS[q]


def within_projector(P, pi):
    C = len(pi); R = np.zeros((C * C, C * C)); V = np.zeros((C * C, C))
    for a in range(C):
        s = np.sqrt(P[a]); R[a * C:(a + 1) * C, a * C:(a + 1) * C] = np.eye(C) - np.outer(s, s)
        for b in range(C): V[a * C + b, b] = math.sqrt(pi[a] * P[a, b])
    return R - R @ V @ np.linalg.pinv(V.T @ R @ V) @ V.T @ R


def basis(G, P, pi):
    """Loewdin-orthonormal basis of the within space from character pairs; returns B (C^2, D) and dir->group index"""
    C = G.C; PW = within_projector(P, pi); w = np.sqrt(pi[:, None] * P).ravel()
    vecs, dg = [], []
    for gi, g in enumerate(G.groups):
        F = np.outer(G.tab[g['i']], np.conj(G.tab[g['j']]))
        parts = [F.real] if g['n'] == 1 else [F.real, F.imag]
        for Fp in parts: vecs.append(PW @ (Fp.ravel() * w)); dg.append(gi)
    Vm = np.array(vecs).T
    Gm = Vm.T @ Vm; e, U = np.linalg.eigh(Gm); S = U @ np.diag(e ** -0.5) @ U.T
    B = Vm @ S
    assert B.shape[1] == (C - 1) ** 2 and np.allclose(B.T @ B, np.eye(B.shape[1]), atol=1e-8)
    assert np.allclose(B @ B.T, PW, atol=1e-6)
    return B, np.array(dg)


def yvec(P, pi):
    C = len(pi); Y = np.zeros((C, C, C * C))
    for a in range(C):
        for b in range(C):
            for b2 in range(C):
                Y[a, b, a * C + b2] = ((1.0 if b2 == b else 0.0) - P[a, b2]) / math.sqrt(pi[a] * P[a, b2])
    return Y


def parse_cells(spec):
    if not spec: return [(q, k) for q in QS_ALL for k in (1, 2)]
    out = []
    for t in spec.split(','):
        q, k = t.split(':'); out.append((int(q), int(k)))
    return out


# ------------------------------------------------------------------ measure
def block_M(X):
    """X (NB, C, C) block tables -> (M = sum_b y y^T, T)"""
    NB, C, _ = X.shape
    T = X.sum(0); P = T / T.sum(1, keepdims=True)
    Na = X.sum(2); E = Na[:, :, None] * P[None]
    with np.errstate(invalid='ignore', divide='ignore'):
        y = np.where(E > 0, (X - E) / np.sqrt(np.where(E > 0, E, 1.0)), 0.0).reshape(NB, C * C)
    return y.T @ y, T


def measure_window(x0, cells, control, rng):
    t0 = time.time()
    P = primes_around(x0)
    pos = P - x0
    res = {}
    bperm = rng.permutation(SPAN // HS[0]) if control == 'bshuf' else np.arange(SPAN // HS[0])
    for q in sorted(set(q for q, _ in cells)):
        G = group(q); C = G.C; cls = G.lut[P % q]
        for k in sorted(k for qq, k in cells if qq == q):
            a = cls[:-k]; b = cls[k:].copy(); ps = pos[:-k]
            inw = (ps >= 0) & (ps < SPAN)
            assert (a[inw] >= 0).all() and (b[inw] >= 0).all()
            aw = a[inw]; bw = b[inw]; pw = ps[inw]
            if control == 'perm':      # successors permuted among pairs with the same first class, whole window
                for c in range(C):
                    ix = np.flatnonzero(aw == c); bw[ix] = bw[ix][rng.permutation(len(ix))]
            ab = aw * C + bw
            for H in HS:
                NB = SPAN // H
                blk = bperm[pw // HS[0]] // (H // HS[0])   # control 'bshuf': the 1e4-blocks in random order
                X = np.bincount(blk * (C * C) + ab, minlength=NB * C * C).reshape(NB, C, C).astype(float)
                M, T = block_M(X)
                res['M_%d_%d_%d' % (q, k, H)] = M; res['NB_%d_%d_%d' % (q, k, H)] = np.array([NB])
            res['T_%d_%d' % (q, k)] = T.astype(np.int64)
    res['x0'] = np.array([x0]); res['lam'] = np.array([math.log(x0 + SPAN / 2)])
    res['nprimes'] = np.array([int(np.sum((pos >= 0) & (pos < SPAN)))])
    return res, time.time() - t0


def cmd_measure(a):
    os.makedirs(a.outdir, exist_ok=True)
    cells = parse_cells(a.cells)
    rng = np.random.default_rng(a.seed)
    for x in a.x0:
        x0 = int(float(x))
        f = os.path.join(a.outdir, 'm_%d.npz' % x0)
        if os.path.exists(f) and not a.force:
            log('skip %d (exists)' % x0); continue
        res, dt = measure_window(x0, cells, a.control, rng)
        res['control'] = np.array([a.control])
        np.savez_compressed(f, **res)
        log('window %d: %d primes, %.1f s' % (x0, res['nprimes'][0], dt))


# ------------------------------------------------------------------ measurement sets
def parse_set(s):
    name, dirs = s.split('=', 1)
    files = []
    for d in dirs.split(','):
        files += sorted(glob.glob(os.path.join(d, 'm_*.npz')))
    assert files, 'no windows in ' + dirs
    return name, files


class MSet:
    def __init__(self, spec):
        self.name, self.files = parse_set(spec)
        self.data = [dict(np.load(f)) for f in self.files]
        self.x0 = [int(d['x0'][0]) for d in self.data]
        self.lam = float(np.mean([d['lam'][0] for d in self.data]))
        self.nw = len(self.data)
        keys = set(self.data[0])
        for d in self.data: keys &= set(d)
        self.cells = sorted(set((int(k.split('_')[1]), int(k.split('_')[2])) for k in keys if k.startswith('T_')))
        self._b = {}

    def pooled(self, q, k):
        T = sum(d['T_%d_%d' % (q, k)] for d in self.data).astype(float)
        return T / T.sum(1, keepdims=True), T.sum(1) / T.sum()

    def basis(self, q, k):
        if (q, k) not in self._b:
            P, pi = self.pooled(q, k); B, dg = basis(group(q), P, pi)
            self._b[(q, k)] = (B, dg, P, pi)
        return self._b[(q, k)]

    def W(self, q, k):
        """per window, H, direction: W_d(H)  -> array (nw, 3, D)"""
        B = self.basis(q, k)[0]
        out = np.zeros((self.nw, len(HS), B.shape[1]))
        for w, d in enumerate(self.data):
            for hi, H in enumerate(HS):
                M = d['M_%d_%d_%d' % (q, k, H)]; NB = int(d['NB_%d_%d_%d' % (q, k, H)][0])
                out[w, hi] = np.einsum('id,ij,jd->d', B, M, B) / (NB - 1)
        return out


# ------------------------------------------------------------------ pair factors
_PR = None
def S_tail(hmax, y):
    """S_(>y)(h), h = 0..hmax (the review's tailemp.py)"""
    global _PR
    if _PR is None: _PR = small_primes(10 ** 7)
    PR = _PR
    p = PR[PR > y].astype(float); Ct = math.exp(np.sum(np.log1p(-1 / (p - 1) ** 2)) - 1 / (1e7 * math.log(1e7)))
    f = np.ones(hmax + 1)
    for pp in PR[(PR > y) & (PR <= hmax)].tolist(): f[pp::pp] *= (pp - 1) / (pp - 2)
    return f * Ct


# ------------------------------------------------------------------ respond
def rough_index(Q, y):
    ok = np.ones(Q, bool)
    for p in small_primes(y).tolist(): ok[::p] = False
    RU = np.flatnonzero(ok); RIDX = -np.ones(Q, np.int64); RIDX[RU] = np.arange(len(RU))
    return RU, RIDX


def count_stretch(x0, q, ks, parts, acc, CH=5 * 10 ** 6, contiguous=True):
    """accumulate per part: n[u]; k = 1: counts (u, cL, cR); k >= 2: signed pair sums (u, a, b)"""
    G = group(q); C = G.C; Q = math.lcm(510510, q)
    RU, RIDX = rough_index(Q, 17); NR = len(RU)
    P = primes_around(x0)
    kmax = max(ks)
    for ci, a in enumerate(range(x0, x0 + SPAN, CH)):
        part = (ci * CH * parts) // SPAN if contiguous else ci % parts
        xs = a + np.arange(CH, dtype=np.int64)
        u = RIDX[xs % Q]; keep = u >= 0; xs = xs[keep]; u = u[keep]
        j = np.searchsorted(P, xs); isp = P[j] == xs
        assert j.min() >= kmax and j.max() + kmax + 1 < len(P)
        cx = G.lut[xs % q]
        off = j - 1 + isp
        Lcls = lambda i: G.lut[P[j - i] % q]
        Rcls = lambda i: G.lut[P[off + i] % q]
        acc['n'][part] += np.bincount(u, minlength=NR)
        base = u * (C * C)
        for k in ks:
            if k == 1:
                acc[1][part] += np.bincount(base + Lcls(1) * C + Rcls(1), minlength=NR * C * C).astype(acc[1].dtype)
            else:
                s = np.bincount(base + Lcls(k) * C + cx, minlength=NR * C * C)
                s += np.bincount(base + cx * C + Rcls(k), minlength=NR * C * C)
                for i in range(1, k):
                    s += np.bincount(base + Lcls(i) * C + Rcls(k - i), minlength=NR * C * C)
                for i in range(1, k + 1):
                    s -= np.bincount(base + Lcls(i) * C + Rcls(k + 1 - i), minlength=NR * C * C)
                acc[k][part] += s.astype(acc[k].dtype)
    return math.log(x0 + SPAN / 2)


def signed_from_counts(cnt, cxu, C):
    """k = 1 counts (NR, C, C) of (cL, cR) -> signed pair sums S (NR, C, C): +(cL,cx) +(cx,cR) -(cL,cR)"""
    S = -cnt.astype(np.int32)
    ar = np.arange(cnt.shape[0])
    S[ar, :, cxu] += cnt.sum(2)
    S[ar, cxu, :] += cnt.sum(1)
    return S


def lr_project(q, k, acc, lam_resp, mset_basis, lam_meas, parts, ys=(11, 13, 17), rev=True, dbatch=12):
    """linear-response tails per direction for one measurement basis"""
    G = group(q); C = G.C; B, dg, P, pi = mset_basis; D = B.shape[1]
    F = (yvec(P, pi) @ B)                                    # (C, C, D)
    Fm = F.reshape(C * C, D)
    Q17 = math.lcm(510510, q); RU, _ = rough_index(Q17, 17); NR = len(RU); cxu = G.lut[RU % q]
    n = acc['n'].astype(float)                               # (parts, NR)
    if k == 1:
        Sp = [signed_from_counts(acc[1][p].reshape(NR, C, C), cxu, C).reshape(NR, C * C) for p in range(parts)]
    else:
        Sp = [acc[k][p].reshape(NR, C * C) for p in range(parts)]
    scale = lam_resp / lam_meas
    Hmax = max(HS); h = np.arange(1, Hmax)
    out = {'scale': scale, 'D': D}
    for y in ys:
        if maxpf(q) > y: continue      # q has a prime factor above y: the comparator is not defined
        Qy = math.lcm(YQ[y], q); RUy, RIDXy = rough_index(Qy, y); NRy = len(RUy)
        inv = RIDXy[RU % Qy]; assert (inv >= 0).all()
        if Qy == Q17:
            Sy = Sp; ny = n
        else:
            Sy = []
            for p in range(parts):
                agg = np.zeros((NRy, C * C), np.int64); np.add.at(agg, inv, Sp[p]); Sy.append(agg.astype(np.int32))
            ny = np.array([np.bincount(inv, weights=n[p], minlength=NRy) for p in range(parts)])
        nempty = int((ny == 0).sum()); nyd = np.maximum(ny, 1)
        St = S_tail(Hmax, y)[1:Hmax]
        wts = {H: ((1 - h[:H - 1] / H) * (1 - np.exp(-h[:H - 1] / lam_meas)) * (St[:H - 1] - 1)) for H in HS}
        hmod = h % Qy
        pref = 2 * (Qy / NRy) / lam_meas * scale
        Wfull = np.zeros((len(HS), D)); Wjk = np.zeros((parts, len(HS), D)); V = np.zeros(D); Vjk = np.zeros((parts, D))
        Wcr = np.zeros((len(HS), D)); Vcr = np.zeros(D)
        ntot = ny.sum(0); nyt = np.maximum(ntot, 1)
        def tail(m1col, A0):
            full = np.zeros(Qy); full[RUy] = m1col; Ff = np.fft.rfft(full)
            A = np.fft.irfft((Ff * Ff.conj()).real, n=Qy) / NRy; A[0] = A0
            Ah = A[hmod]
            return [pref * np.dot(wts[H], Ah[:H - 1]) for H in HS]
        for d0 in range(0, D, dbatch):
            dd = slice(d0, min(D, d0 + dbatch))
            num = [(Sy[p] @ Fm[:, dd]) for p in range(parts)]                      # (NRy, nd) per part
            m1 = [num[p] / nyd[p][:, None] for p in range(parts)]
            numt = sum(num)
            for jd in range(numt.shape[1]):
                d = d0 + jd
                mp = np.array([m1[p][:, jd] for p in range(parts)])                 # (parts, NRy)
                sm = mp.sum(0); s2 = (mp ** 2).sum(0)
                # lag 0: cross products of disjoint parts (no sampling-noise term)
                A0 = ((sm ** 2 - s2).sum() / (parts * (parts - 1))) / NRy
                Wfull[:, d] = tail(numt[:, jd] / nyt, A0); V[d] = A0 * scale
                for L in range(parts):
                    smL = sm - mp[L]; s2L = s2 - mp[L] ** 2
                    A0L = ((smL ** 2 - s2L).sum() / ((parts - 1) * (parts - 2))) / NRy
                    nL = np.maximum(ntot - ny[L], 1)
                    Wjk[L, :, d] = tail((numt[:, jd] - num[L][:, jd]) / nL, A0L); Vjk[L, d] = A0L * scale
                # all-lag cross-part estimate (systematic check of the estimator)
                Fp = []
                for p in range(parts):
                    full = np.zeros(Qy); full[RUy] = mp[p]; Fp.append(np.fft.rfft(full))
                sF = sum(Fp); sF2 = sum((f * f.conj()).real for f in Fp)
                A = np.fft.irfft(((sF * sF.conj()).real - sF2) / (parts * (parts - 1)), n=Qy) / NRy
                Ah = A[hmod]; Wcr[:, d] = [pref * np.dot(wts[H], Ah[:H - 1]) for H in HS]; Vcr[d] = A[0] * scale
        out['y%d' % y] = dict(W=Wfull.tolist(), Wjk=Wjk.tolist(), V=V.tolist(), Vjk=Vjk.tolist(), Wcross=Wcr.tolist(),
                              Vcross=Vcr.tolist(), Q=Qy, nclasses=NRy,
                              nmin=float(ny.sum(0).min()), nmean=float(ny.sum(0).mean()), nempty_classparts=nempty)
    if rev and k == 1 and C <= 6:     # memory: only for the small moduli (validation against the review)
        # the review's estimate: pooled counts, per-class mean and the noise correction var/(n-1) at lag 0
        cnt = sum(acc[1][p] for p in range(parts)).reshape(NR, C, C).astype(float); nn = cnt.sum((1, 2))
        RQ = np.zeros((C, C, C, D))                          # rq[cx, cL, cR, d]
        for c in range(C): RQ[c] = F[:, c, :][:, None, :] + F[c][None, :, :] - F
        m1 = np.zeros((NR, D)); m2 = np.zeros((NR, D))
        for c in range(C):
            sel = cxu == c; cc = cnt[sel].reshape(-1, C * C)
            m1[sel] = cc @ RQ[c].reshape(C * C, D); m2[sel] = cc @ (RQ[c].reshape(C * C, D) ** 2)
        m1 /= nn[:, None]; m2 /= nn[:, None]
        noise = (m2 - m1 ** 2) / (nn[:, None] - 1)
        Vrev = ((m1 ** 2) - noise).mean(0) * scale
        St = S_tail(Hmax, 17)[1:Hmax]; hmod = h % Q17; pref = 2 * (Q17 / NR) / lam_meas * scale
        Wrev = np.zeros((len(HS), D))
        for d in range(D):
            full = np.zeros(Q17); full[RU] = m1[:, d]; Ff = np.fft.rfft(full)
            A = np.fft.irfft((Ff * Ff.conj()).real, n=Q17) / NR; A[0] -= noise[:, d].mean()
            Ah = A[hmod]
            for hi, H in enumerate(HS):
                Wrev[hi, d] = pref * np.dot((1 - h[:H - 1] / H) * (1 - np.exp(-h[:H - 1] / lam_meas)) * (St[:H - 1] - 1), Ah[:H - 1])
        out['rev'] = dict(V=Vrev.tolist(), W=Wrev.tolist())
    return out


def cmd_respond(a):
    sets = [MSet(s) for s in a.set]
    ks = [int(v) for v in a.ks.split(',')]
    parts = a.parts
    out = json.load(open(a.out)) if (os.path.exists(a.out) and a.append) else {}
    for q in [int(v) for v in a.qs.split(',')]:
        G = group(q); C = G.C; Q = math.lcm(510510, q); RU, _ = rough_index(Q, 17); NR = len(RU)
        acc = {'n': np.zeros((parts, NR), np.int64)}
        for k in ks: acc[k] = np.zeros((parts, NR * C * C), np.int32)
        lams = []; stretches = []
        for f in (a.raw_in or []):
            R = np.load(f)
            if 'n_%d' % q not in R: continue
            acc['n'] += R['n_%d' % q]
            for k in ks: acc[k] += R['acc_%d_%d' % (q, k)]
            lams += list(R['lams_%d' % q]); stretches += list(R['stretches_%d' % q])
        for x in (a.stretch or []):
            x0 = int(float(x)); t0 = time.time()
            lams.append(count_stretch(x0, q, ks, parts, acc)); stretches.append(x0)
            log('q=%d: counted stretch %d in %.1f s' % (q, x0, time.time() - t0))
        if a.raw_out:
            fn = a.raw_out.replace('.npz', '_q%d.npz' % q)
            np.savez_compressed(fn, **{'n_%d' % q: acc['n'], 'lams_%d' % q: np.array(lams),
                                       'stretches_%d' % q: np.array(stretches), **{'acc_%d_%d' % (q, k): acc[k] for k in ks}})
            log('raw counts written to %s' % fn)
        lam_resp = float(np.mean(lams))
        for S in sets:
            for k in ks:
                if (q, k) not in S.cells: continue
                t0 = time.time()
                r = lr_project(q, k, acc, lam_resp, S.basis(q, k), S.lam, parts, rev=not a.no_rev)
                r.update(lam_resp=lam_resp, lam_meas=S.lam, stretches=[int(v) for v in stretches], set=S.name,
                         nwin=S.nw, P=S.basis(q, k)[2].tolist())
                out['%s|%d_%d' % (S.name, q, k)] = r
                log('q=%d k=%d set %s: projected in %.1f s' % (q, k, S.name, time.time() - t0))
                json.dump(out, open(a.out, 'w'))
        del acc


# ------------------------------------------------------------------ chain null
def cmd_chain(a):
    S = MSet(a.set[0])
    cells = parse_cells(a.cells)
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    lam = S.lam; nw = a.nwin; NB4 = SPAN // HS[0]; per = HS[-1] // HS[0]      # 100 sub-blocks per 1e6 block
    nchain = nw * (SPAN // HS[-1])
    res = [dict() for _ in range(nw)]
    for q in sorted(set(q for q, _ in cells)):
        G = group(q); C = G.C; ks = sorted(k for qq, k in cells if qq == q); kmax = max(ks)
        P, _ = S.pooled(q, 1)
        w, v = np.linalg.eig(P.T); pi = np.real(v[:, np.argmin(abs(w - 1))]); pi /= pi.sum()
        cum = np.cumsum(P, 1); cum[:, -1] = 1.0
        sub = rng.poisson(HS[0] / lam, (nchain, per))                         # transitions per 1e4 sub-block
        nb = sub.sum(1); L = int(nb.max()) + kmax + 1
        t0 = time.time()
        s = np.empty((nchain, L), np.int8)
        s[:, 0] = rng.choice(C, nchain, p=pi)
        for t in range(1, L):
            uu = rng.random(nchain)
            s[:, t] = (uu[:, None] > cum[s[:, t - 1].astype(np.int64)]).sum(1)
        log('q=%d: simulated %d chains of %d steps in %.1f s' % (q, nchain, L, time.time() - t0))
        csub = np.cumsum(sub, 1)
        for wi in range(nw):
            rows = range(wi * (SPAN // HS[-1]), (wi + 1) * (SPAN // HS[-1]))
            for k in ks:
                ab = []; b4 = []
                for bi, r in enumerate(rows):
                    n_ = nb[r]; t = np.arange(n_)
                    ab.append(s[r, :n_].astype(np.int64) * C + s[r, k:k + n_])
                    b4.append(bi * per + np.searchsorted(csub[r], t, side='right'))
                ab = np.concatenate(ab); b4 = np.concatenate(b4)
                for H in HS:
                    NB = SPAN // H; blk = b4 // (H // HS[0])
                    X = np.bincount(blk * (C * C) + ab, minlength=NB * C * C).reshape(NB, C, C).astype(float)
                    M, T = block_M(X)
                    res[wi]['M_%d_%d_%d' % (q, k, H)] = M; res[wi]['NB_%d_%d_%d' % (q, k, H)] = np.array([NB])
                res[wi]['T_%d_%d' % (q, k)] = T.astype(np.int64)
        del s
    for wi in range(nw):
        idx = a.first + wi
        res[wi]['x0'] = np.array([-(idx + 1)]); res[wi]['lam'] = np.array([lam]); res[wi]['control'] = np.array(['chain'])
        np.savez_compressed(os.path.join(a.outdir, 'm_%d.npz' % idx), **res[wi])
    log('chain windows %d..%d written' % (a.first, a.first + nw - 1))


# ------------------------------------------------------------------ report
def jk_se(vals):
    """jackknife standard error of the mean over windows (rows)"""
    v = np.asarray(vals, float); n = len(v)
    loo = (v.sum(0) - v) / (n - 1)
    return np.sqrt((n - 1) / n * ((loo - loo.mean(0)) ** 2).sum(0))


def summarise(S, q, k, lr):
    G = group(q); B, dg, P, pi = S.basis(q, k); D = B.shape[1]
    Wd = S.W(q, k)                                            # (nw, 3, D)
    sl = (Wd[:, 2, :] - Wd[:, 0, :]) / LN100                  # (nw, D)
    have = lr is not None
    if have:
        assert np.allclose(np.array(lr['P']), P), 'LR file was projected on another basis'
        y = lr['y17']; Wt = np.array(y['W']); Wjk = np.array(y['Wjk'])
        lrs = (Wt[2] - Wt[0]) / LN100; lrs_jk = (Wjk[:, 2] - Wjk[:, 0]) / LN100
        asym = -np.array(y['V']) / S.lam
        def lrs_of(key, W='W'):
            if key not in lr or maxpf(q) > int(key[1:]): return np.full(D, np.nan)
            Wy = np.array(lr[key][W]); return (Wy[2] - Wy[0]) / LN100
        lr13 = lrs_of('y13'); lr11 = lrs_of('y11'); lrcross = lrs_of('y17', 'Wcross')
        rev = lr.get('rev')
        revasym = -np.array(rev['V']) / S.lam if rev else None
    def agg(sel, name, extra):
        sel = np.asarray(sel)
        Wm = Wd[:, :, sel].mean(2)                            # (nw, 3)
        s = sl[:, sel].mean(1)
        r = dict(name=name, ndir=int(len(sel)), **extra,
                 W=[[float(Wm[:, i].mean()), float(Wm[:, i].std(ddof=1) / math.sqrt(S.nw))] for i in range(3)],
                 slope=[float(s.mean()), float(jk_se(s))], slope_pw=s.tolist())
        if have:
            l = float(lrs[sel].mean()); ljk = lrs_jk[:, sel].mean(1); np_ = len(ljk)
            le = float(math.sqrt((np_ - 1) / np_ * ((ljk - ljk.mean()) ** 2).sum()))
            em = float(jk_se(s - l))
            r.update(lr=l, lr_err=le, lr_asym=float(asym[sel].mean()), lr13=float(lr13[sel].mean()), lr11=float(lr11[sel].mean()),
                     lr_cross=float(lrcross[sel].mean()),
                     excess=float(s.mean() - l), err_meas=em, err_tot=float(math.hypot(em, le)),
                     z=float((s.mean() - l) / math.hypot(em, le)),
                     excess_asym=float(s.mean() - asym[sel].mean()),
                     V=float(np.array(y['V'])[sel].mean()))
            if revasym is not None: r.update(lr_asym_rev=float(revasym[sel].mean()), V_rev=float(np.array(rev['V'])[sel].mean()))
        return r
    groups = []
    for gi, g in enumerate(G.groups):
        groups.append(agg(np.flatnonzero(dg == gi), g['label'], dict(type=g['type'], sub=g['sub'], cond=g['cond'])))
    types = {}
    for t in ('eq', 'conj', 'rr', 'rc', 'cc'):
        sel = [d for d in range(D) if G.groups[dg[d]]['type'] == t]
        if sel: types[t] = agg(sel, t, {})
    subs = {}
    for key in sorted(set((g['type'], g['sub']) for g in G.groups)):
        sel = [d for d in range(D) if (G.groups[dg[d]]['type'], G.groups[dg[d]]['sub']) == key]
        subs['%s %s' % key] = agg(sel, '%s %s' % key, {})
    tot = agg(list(range(D)), 'all', {})
    return dict(q=q, k=k, D=D, lam=S.lam, nw=S.nw, groups=groups, types=types, subtypes=subs, all=tot, has_lr=have)


def fmt_row(r, have):
    s = '  %-22s %3d  W %.4f(%2.0f) %.4f(%2.0f) %.4f(%3.0f)  slope %+.4f(%2.0f)' % (
        r['name'], r['ndir'], r['W'][0][0], r['W'][0][1] * 1e4, r['W'][1][0], r['W'][1][1] * 1e4, r['W'][2][0], r['W'][2][1] * 1e4,
        r['slope'][0], r['slope'][1] * 1e4)
    if have and 'lr' in r:
        s += ' | LR %+.4f(%2.0f) asym %+.4f | excess %+.4f +- %.4f (%+.1f)  [y13 %+.4f y11 %+.4f]' % (
            r['lr'], r['lr_err'] * 1e4, r['lr_asym'], r['excess'], r['err_tot'], r['z'], r['lr13'], r['lr11'])
    return s


def cmd_report(a):
    S = MSet(a.set[0])
    LR = {}
    for f in (a.lr or []):
        for key, v in json.load(open(f)).items():
            if key.split('|')[0] == S.name: LR[key.split('|')[1]] = v
    out = dict(set=S.name, windows=S.x0, lam=S.lam, cells={})
    log('set %s: %d windows, lam %.4f   (numbers in brackets: error in units of 1e-4)' % (S.name, S.nw, S.lam))
    for q, k in S.cells:
        if a.qs and q not in [int(v) for v in a.qs.split(',')]: continue
        lr = LR.get('%d_%d' % (q, k))
        r = summarise(S, q, k, lr); out['cells']['%d_%d' % (q, k)] = r
        log('q=%d k=%d  D=%d  %s' % (q, k, r['D'], '' if lr else '(no LR)'))
        log(fmt_row(r['all'], r['has_lr']))
        for t, v in r['types'].items(): log(fmt_row(v, r['has_lr']))
        if a.groups:
            for v in r['subtypes'].values(): log('   sub' + fmt_row(v, r['has_lr']))
            for v in r['groups']: log('   grp' + fmt_row(v, r['has_lr']) + '  [%s %s cond %s]' % (v['type'], v['sub'], v['cond']))
    if a.json: json.dump(out, open(a.json, 'w'), indent=0)


def cmd_combine(a):
    R = [json.load(open(f)) for f in a.reports]
    log('combine %s' % ' + '.join('%s (lam %.2f, %d windows)' % (r['set'], r['lam'], len(r['windows'])) for r in R))
    out = {}
    for cell in R[0]['cells']:
        if not all(cell in r['cells'] for r in R): continue
        rows = {}
        for kind in ('types', 'subtypes'):
            for name in R[0]['cells'][cell][kind]:
                xs = [r['cells'][cell][kind][name] for r in R]
                if not all('excess' in x for x in xs): continue
                w = np.array([1 / x['err_tot'] ** 2 for x in xs]); e = np.array([x['excess'] for x in xs])
                m = float((w * e).sum() / w.sum()); se = float(1 / math.sqrt(w.sum()))
                rows[name] = dict(excess=m, err=se, z=m / se, parts=[[x['excess'], x['err_tot']] for x in xs])
                log('  %s %-14s excess %+.4f +- %.4f (%+.1f)   parts %s' % (cell, name, m, se, m / se,
                    '  '.join('%+.4f+-%.4f' % (x['excess'], x['err_tot']) for x in xs)))
        out[cell] = rows
    if a.json: json.dump(out, open(a.json, 'w'), indent=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    m = sp.add_parser('measure'); m.add_argument('outdir'); m.add_argument('x0', nargs='+')
    m.add_argument('--cells', default=''); m.add_argument('--control', default='none', choices=['none', 'perm', 'bshuf'])
    m.add_argument('--seed', type=int, default=1); m.add_argument('--force', action='store_true')
    r = sp.add_parser('respond'); r.add_argument('--out', required=True); r.add_argument('--stretch', nargs='*')
    r.add_argument('--set', action='append', default=[]); r.add_argument('--qs', default='3,5,7,8,9,11,12,13')
    r.add_argument('--ks', default='1,2'); r.add_argument('--parts', type=int, default=4)
    r.add_argument('--raw-out'); r.add_argument('--raw-in', nargs='*'); r.add_argument('--append', action='store_true')
    r.add_argument('--no-rev', action='store_true')
    c = sp.add_parser('chain'); c.add_argument('outdir'); c.add_argument('--set', action='append', required=True)
    c.add_argument('--nwin', type=int, default=24); c.add_argument('--seed', type=int, default=1)
    c.add_argument('--cells', default=''); c.add_argument('--first', type=int, default=0)
    p = sp.add_parser('report'); p.add_argument('--set', action='append', required=True); p.add_argument('--lr', nargs='*')
    p.add_argument('--json'); p.add_argument('--groups', action='store_true'); p.add_argument('--qs', default='')
    cb = sp.add_parser('combine'); cb.add_argument('reports', nargs='+'); cb.add_argument('--json')
    a = ap.parse_args()
    dict(measure=cmd_measure, respond=cmd_respond, chain=cmd_chain, report=cmd_report, combine=cmd_combine)[a.cmd](a)


if __name__ == '__main__':
    main()
