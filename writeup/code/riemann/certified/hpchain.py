"""High-precision certified chain: the decoder of lpchain.py on samples computed in ball arithmetic (hpwin.py),
with every proof replayed in ball arithmetic before a bit is committed.

- Samples: zero side from Arb's rigorous zeros, kernel and transform as 0F1 balls, trivial zeros integrated rigorously;
  the bound delta covers only the zeros above T (no precision or rounding terms are needed: the balls carry them).
- Search: lpchain's propagation, Farkas LPs and branch and bound, in double precision on the ball midpoints
  (bounds rounded outward). Search only proposes; nothing is trusted from it.
- Proof: each certified bit's tree (forcings, violated rows, Farkas vectors, splits) is replayed with the balls:
  a forcing or violation must hold for every point of every ball, a Farkas vector must satisfy
  lower(min over the box of (A^T w).b) > upper(w+.hi - w-.lo). A bit whose replay fails is not committed.
- easy_first: certify bits anywhere in the window with growing node budgets, keep every proved bit, repeat.
Usage: hpchain.py X stages opts   (as lpchain.py; opts may set easy_first, budgets, maxpos)
"""
import sys, os, math, json, time
import numpy as np
from flint import arb, arb_mat
import lpchain
from lpchain import propagate, farkas, refute
import hpwin

lpchain.REL = float(os.environ.get('REL', '1e-15'))   # the search may propose thin certificates; the replay decides

def down(v): return float(np.nextafter(float(v), -np.inf))
def up(v): return float(np.nextafter(float(v), np.inf))

class HPChain(lpchain.LPChain):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.zcache = {}; self.replays = dict(ok=0, failed=0); self.white = True

    def zs(self, x, a, beta, m, c):
        key = (x, a, beta, m, c)
        v = self.zcache.get(key)
        if v is None:
            T = c * beta * x / a
            self.maxT = max(self.maxT, T)
            v = self.zcache[key] = hpwin.sample(x, a, beta, m, T)
        return v

    def kern(self, x, n, a, beta, m):
        xa = arb(x); v = xa * (arb(n) / xa).log()
        w = 1 - (v / arb(a)) ** 2
        if w.upper() <= 0: return arb(0)
        if w.lower() <= 0:   # touches the edge: K is increasing in w and 0 there
            hi = hpwin.K(arb(a) * (1 - w.upper()).sqrt(), a, beta, m).upper()
            return arb(0).union(arb(hi))
        return hpwin.K(v, a, beta, m)

    def window(self, c, a, beta, m, cc, k, h):
        U, w = self.unknowns(c, k)
        xs = []
        x = max(2.0, math.floor(c - a - 1))
        while x < w:
            if x * math.exp(-a / x) > 1.0 and abs(x * math.log(c / x)) < a and x * math.log(w / x) >= a:
                xs.append(float(x))
            x += h
        k0 = hpwin.K0(beta, m)
        Uset = set(U)
        Aa = []; loa = []; hia = []
        for x in xs:
            z, d = self.zs(x, a, beta, m, cc)
            lo_, hi_ = math.floor(x * math.exp(-a / x)) - 1, math.ceil(x * math.exp(a / x)) + 1
            known = arb(0)
            for n in range(max(2, lo_), hi_ + 1):
                if n in Uset: continue
                lv = self.lam[n] if n < len(self.lam) else np.nan
                kv = self.kern(x, n, a, beta, m)
                if kv.upper() == 0: continue
                if np.isnan(lv):
                    raise RuntimeError(f'unknown number {n} inside the support of sample {x}')
                if lv > 0:
                    p = round(math.exp(lv))
                    known += arb(p).log() * (arb(x) / n).sqrt() * kv
            y = (z - known) / k0; dd = arb(d) / k0
            loa.append((y - dd).lower()); hia.append((y + dd).upper())
            Aa.append([arb(u).log() * (arb(x) / u).sqrt() * self.kern(x, u, a, beta, m) / k0 for u in U])
        n0 = len(Aa); self.n0 = n0
        A0 = np.array([[float(v.mid()) for v in row] for row in Aa])
        self.weights = []; self.scales = []
        AM0 = arb_mat(Aa)
        if self.white:
            # implied rows u^T (A b) in [u+ . lo - u- . hi, u+ . hi - u- . lo], u = left singular vectors (double, used as exact weights)
            Uw, sv, _ = np.linalg.svd(A0, full_matrices=False)
            r = int(np.sum(sv >= 1e-15 * sv[0]))
            if r:
                W = Uw[:, :r].T
                Wp = np.where(W > 0, W, 0.0); Wm = np.where(W < 0, W, 0.0)
                rowsM = arb_mat(W.tolist()) * AM0
                loM = arb_mat([[v] for v in loa]); hiM = arb_mat([[v] for v in hia])
                lw = arb_mat(Wp.tolist()) * loM + arb_mat(Wm.tolist()) * hiM
                hw = arb_mat(Wp.tolist()) * hiM + arb_mat(Wm.tolist()) * loM
                for i in range(r):
                    row = [rowsM[i, j] for j in range(len(U))]
                    scf = 1.0 / max(abs(float(v.mid())) for v in row); sc = arb(scf)
                    self.weights.append(W[i].tolist()); self.scales.append(scf)
                    Aa.append([v * sc for v in row]); loa.append((lw[i, 0] * sc).lower()); hia.append((hw[i, 0] * sc).upper())
        self.last_xs = xs; self.Aa = Aa; self.loa = loa; self.hia = hia
        A = np.array([[float(v.mid()) for v in row] for row in Aa])
        lo = np.array([down(v) for v in loa]); hi = np.array([up(v) for v in hia])
        self.bits_info = float((arb(math.log(c)) / ((hia[len(hia) // 2] - loa[len(loa) // 2]) / 2)).log().mid() / math.log(2))
        return U, A, lo, hi

    # ---- replay of one refutation tree in ball arithmetic ----
    def replay(self, nodes, kk):
        """Replay a refutation tree with the balls (Arb matrices): every event must hold for every point of every ball."""
        AM = arb_mat(self.Aa); lo, hi = self.loa, self.hia; n = len(self.Aa)
        dec = lambda s: [{'x': -1, '0': 0, '1': 1}[ch] for ch in s]
        def sums(st):
            e1 = arb_mat([[1 if v == 1 else 0] for v in st]); ef = arb_mat([[1 if v == -1 else 0] for v in st])
            return AM * e1, AM * ef
        states = set(); need = []
        for node in nodes:
            st = dec(node['state']); states.add(node['state']); closed = False
            for ev in node['events']:
                base, tot = sums(st)
                if ev[0] == 'force':
                    _, var, v, row = ev
                    ok = (base[row, 0] + AM[row, var]).lower() > hi[row] if v == 0 else (base[row, 0] + tot[row, 0] - AM[row, var]).upper() < lo[row]
                    if st[var] != -1 or not ok: return False
                    st[var] = v
                else:
                    _, row, side = ev
                    ok = base[row, 0].lower() > hi[row] if side == 'hi' else (base[row, 0] + tot[row, 0]).upper() < lo[row]
                    if not ok: return False
                    closed = True
            if closed: continue
            if node.get('farkas'):
                w = node['farkas']; free = [q for q in range(kk) if st[q] == -1]
                base, _ = sums(st)
                g = arb_mat([w]) * AM
                minbox = arb(0)
                for q in free:
                    minbox += min(g[0, q].lower(), 0)
                rhs = arb(0)
                for i in range(n):
                    if w[i] == 0: continue
                    rhs += arb(w[i]) * ((hi[i] if w[i] > 0 else lo[i]) - base[i, 0])
                if not minbox.lower() > rhs.upper(): return False
            elif 'branch' in node:
                q = node['branch']
                for v in (0, 1):
                    s2 = list(st); s2[q] = v; need.append(''.join('x01'[t + 1] for t in s2))
            else:
                return False
        return all(s in states for s in need)

    def refute(self, A, lo, hi, st0, limit, log=None):
        n0 = self.n0; A0, lo0, hi0 = A[:n0], lo[:n0], hi[:n0]
        stack = [st0.copy()]; nodes = 0; lps = 0
        while stack:
            st = stack.pop(); nodes += 1
            if nodes > limit: return ('limit', nodes, lps)
            rec = None if log is None else []
            node = None if log is None else dict(state=''.join('x01'[v + 1] for v in st), events=rec)
            if node is not None: log.append(node)
            if not propagate(A0, lo0, hi0, st, rec): continue
            free = np.where(st < 0)[0]
            if len(free) == 0:
                # a full pattern that fits the original rows: check the implied rows too
                tot = A[:, st == 1].sum(1)
                if np.all(tot >= lo - 1e-9 * (1 + np.abs(lo))) and np.all(tot <= hi + 1e-9 * (1 + np.abs(hi))):
                    return ('feasible', st.copy())
                if node is not None: node['farkas'] = None
                return ('limit', nodes, lps)
            base = A[:, st == 1].sum(1)
            lps += 1
            ok, v = farkas(A[:, free], lo - base, hi - base)
            if ok:
                if node is not None:
                    w = np.zeros(A.shape[0]); w[np.where(farkas.keep)[0]] = v; node['farkas'] = w.tolist()
                continue
            j = 0; r = int(round(v[j]))
            if node is not None: node['branch'] = int(free[j])
            for val in (1 - r, r):
                s2 = st.copy(); s2[free[j]] = val; stack.append(s2)
        return ('refuted', nodes, lps)

    def certify(self, A, lo, hi, st, j, val, limit):
        """Refute b_j = 1 - val under the fixed bits st; True only if the double search refutes it and the replay holds."""
        s2 = st.copy(); s2[j] = 1 - val
        lg = []
        rr = self.refute(A, lo, hi, s2, limit, lg)
        if rr[0] != 'feasible':
            self.stats['nodes'] += rr[1]; self.stats['lps'] += rr[2]; self.stats['maxnodes'] = max(self.stats['maxnodes'], rr[1])
        if rr[0] != 'refuted': return rr[0], rr
        ok = self.replay(lg, len(st))
        if ok and self.rec_now is not None:
            self.rec_now.append(dict(bit=int(j), value=int(val), root=''.join('x01'[v + 1] for v in s2), nodes=lg))
        self.replays['ok' if ok else 'failed'] += 1
        return ('refuted' if ok else 'replay failed'), rr

    def decide(self, U, A, lo, hi, limit, maxbits):
        k = len(U); st = np.full(k, -1)
        self.rec_now = [] if self.record and self.record[0] <= U[0] < self.record[1] else None
        if self.rec_now is not None:
            self.proofs[U[0]] = dict(U=[int(u) for u in U], xs=self.last_xs, params=self.cur, n0=self.n0, weights=self.weights, scales=self.scales, proofs=self.rec_now)
        budgets = self.budgets if self.easy else [limit]
        bi = 0
        while bi < len(budgets):
            B = budgets[bi]; progress = False
            s0 = st.copy()
            if not propagate(A[:self.n0], lo[:self.n0], hi[:self.n0], s0):
                return [], 'data inconsistent'
            free = np.where(s0 < 0)[0]
            if len(free):
                base = A[:, s0 == 1].sum(1)
                ok, v = farkas(A[:, free], lo - base, hi - base)
                if ok: return [], 'data inconsistent'
                guess = np.zeros(k, dtype=int); guess[free] = np.rint(v).astype(int); guess[s0 >= 0] = s0[s0 >= 0]
            else:
                guess = s0
            positions = range(min(k, self.maxpos)) if self.easy else [int(np.argmax(st < 0))]
            fits = {}
            for j in positions:
                if st[j] >= 0: continue
                for val in (int(guess[j]), 1 - int(guess[j])):
                    res, rr = self.certify(A, lo, hi, st, j, val, B)
                    if res == 'refuted':
                        st[j] = val; progress = True; break
                    if res == 'feasible': fits[1 - val] = True
                if not self.easy and st[j] < 0:
                    return self.prefix(st), ('ambiguous' if len(fits) == 2 else f'search limit at {U[j]}')
                if progress and B > budgets[0]: break
            lead = self.prefix(st)
            if len(lead) >= maxbits or len(lead) == k: break
            bi = 0 if progress else bi + 1
        got = self.prefix(st)
        self.extra = {U[j]: int(st[j]) for j in range(len(got), k) if st[j] >= 0}
        return got, (None if got else 'search limit')

    def prefix(self, st):
        got = []
        for v in st:
            if v < 0: break
            got.append(int(v))
        return got

    def run(self, stages, start=3, limit=1500, maxbits=40, h=1.0, verbose=False, record=(), easy_first=False, budgets=(10, 40, 160, 640, 2560), maxpos=60, white=True):
        self.easy = easy_first; self.white = white; self.budgets = list(budgets); self.maxpos = maxpos; self.record = tuple(record)
        c = start; t0 = time.time(); decided = 0; self.extra = {}
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
                    return dict(stop=c, reason=why, decided=decided, secs=round(time.time() - t0, 1), bits_at_stop=round(self.bits_info, 1), replays=self.replays, **self.stats)
                for u, v in zip(U, got):
                    self.commit(u, v == 1); decided += 1
                for u, v in self.extra.items():      # proved bits beyond the prefix (easy-first)
                    self.commit(u, v == 1); decided += 1
                self.extra = {}
                c = U[len(got) - 1] + 1
                if verbose:
                    print(c, decided, round(time.time() - t0), f'bits {self.bits_info:.1f}', self.replays, self.stats, flush=True)
            self.log.append(dict(end=end, a=a, beta=beta, m=m, c=cc, k=k, height_over_pi_x=cc * beta / a / math.pi))
        return dict(stop=None, decided=decided, maxT=self.maxT, secs=round(time.time() - t0, 1), replays=self.replays, **self.stats)

if __name__ == '__main__':
    X = int(sys.argv[1]); stages = json.loads(sys.argv[2])
    opts = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    wheel = tuple(opts.pop('wheel', [2])); tb = opts.pop('truth_below', 0); start = opts.pop('start', 3)
    out = opts.pop('out', None)
    C = HPChain(X, wheel, tb)
    r = C.run(stages, start=start, **opts)
    print(r)
    from sympy import isprime
    lo = max(start, 3)
    primes = [n for n in range(lo, X) if C.lam[n] == math.log(n)]
    hi_dec = r['stop'] or X
    wrong = [n for n in range(lo, hi_dec) if not np.isnan(C.lam[n]) and (C.lam[n] == math.log(n)) != isprime(n)]
    print(f'{len(primes)} primes certified in [{lo}, {hi_dec}); wrong: {wrong[:10]}; zeros up to {C.maxT:.1f}')
    if out:
        json.dump(dict(X=X, stages=stages, start=start, result=r, primes=primes, opts=opts), open(out, 'w'))
        if C.proofs:
            json.dump({str(k): v for k, v in C.proofs.items()}, open(out.replace('.json', '_proofs.json'), 'w'))
