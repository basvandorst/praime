"""Round 16 (preregister/round16.txt): nearest neighbours of Eisenstein primes from Hardy-Littlewood in Z[w], w = e^(2 pi i/3).
Adapted from r11.py (Gaussian primes): norm x^2 - xy + y^2, residue maps of Z[w], density 3 sqrt(3)/(pi ln N).
usage: r11.py data | theory | report"""
import sys, os, math, json, time
from itertools import combinations
import numpy as np
import gmpy2
from scipy.stats import t as student_t
from common import small_primes
import theory2 as T2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
BOXES = {'E1': (200_000, 100_000, 2000), 'E2': (6_000_000, 3_000_000, 1000)}
SUF = ''
MARGIN = 40
RMAX = 75
RTH = 52
KAPPA = 3 * math.sqrt(3) / math.pi
NRM = lambda x, y: x * x - x * y + y * y
NORMS = sorted({NRM(x, y) for x in range(-12, 13) for y in range(-12, 13) if 0 < NRM(x, y) <= RTH})
REPORT_D = [1, 3, 4, 7, 9, 12, 13]

_P = small_primes(20_000_000)
_SPLIT = _P[_P % 3 == 1].astype(np.float64)
_INERT = _P[_P % 3 == 2].astype(np.float64)
_INERT = _INERT[_INERT > 2] if False else _INERT


def G_tail(k, P0):
    s = _SPLIT[_SPLIT > P0]; i = _INERT[_INERT > P0]
    i = i[i * i <= 2e7]
    v = 2 * np.sum(np.log1p(-k / s) - k * np.log1p(-1 / s)) + np.sum(np.log1p(-k / i ** 2) - k * np.log1p(-1 / i ** 2))
    return math.exp(v)


_I0 = {}


def i0(p):
    if p not in _I0:
        _I0[p] = next(x for x in range(2, p) if (x * x + x + 1) % p == 0)
    return _I0[p]


def distinct(cls):
    cls = np.sort(cls, axis=1)
    return 1 + np.count_nonzero(np.diff(cls, axis=1), axis=1)


def S_G(sets, P0):
    """sets: int array (n, k, 2). Singular series of Hardy-Littlewood in Z[w]."""
    n, k, _ = sets.shape
    x, y = sets[:, :, 0], sets[:, :, 1]
    logS = np.zeros(n); dead = np.zeros(n, dtype=bool)

    def fac(nu, N):
        f = 1 - nu / N
        nonlocal dead
        dead |= f <= 0
        return np.log(np.where(f > 0, f, 1.0)) - k * math.log(1 - 1 / N)
    logS += fac(distinct((x + y) % 3), 3.0)                               # 1 - w (norm 3): w -> 1
    for p in _P[(_P != 3) & (_P <= P0)]:
        p = int(p)
        if p % 3 == 1:
            r = i0(p); r2 = (-1 - r) % p                                  # the two roots of x^2 + x + 1
            logS += fac(distinct((x + y * r) % p), float(p)) + fac(distinct((x + y * r2) % p), float(p))
        else:
            logS += fac(distinct((x % p) * p + (y % p)), float(p * p))
    S = np.exp(logS) * G_tail(k, P0)
    S[dead] = 0.0
    return S


class Disc:
    """Bahadur pattern: B = {0}, interior = lattice points w with 0 < |w|^2 < R (x + y even)."""

    def __init__(self, R, order):
        pts = [(x, y) for x in range(-12, 13) for y in range(-12, 13) if 0 < NRM(x, y) < R]
        self.SB = 1.0
        I = np.array(pts, dtype=np.int64)
        P0 = max(4 * R, 50)
        zero = np.zeros((1, 1, 2), dtype=np.int64)
        SB = float(S_G(zero, P0)[0])
        Ssub = lambda cols: S_G(np.concatenate([np.repeat(zero, len(cols), 0), cols], axis=1), P0) / SB
        self.s1 = np.zeros(0); self.pairs = self.trip = self.quad = None
        if len(I) == 0:
            return
        s1 = Ssub(I[:, None, :]); live = I[s1 > 0]; s1 = s1[s1 > 0]; n = len(live)
        self.s1 = s1; self.n = n
        if n >= 2:
            i, j = np.triu_indices(n, 1)
            v = Ssub(np.stack([live[i], live[j]], 1))
            self.s2 = np.zeros((n, n)); self.s2[i, j] = v; self.s2[j, i] = v
            self.pairs = (i, j, v)
        if order >= 3 and n >= 3:
            T = np.array(list(combinations(range(n), 3)), dtype=np.int64)
            a, b, c = T.T
            s3v = Ssub(live[T])
            self.trip = (a, b, c, s3v)
            if order >= 4 and n >= 4:
                s3 = np.zeros((n, n, n))
                for (p1, p2, p3) in ((a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a)):
                    s3[p1, p2, p3] = s3v
                Q = np.array(list(combinations(range(n), 4)), dtype=np.int64)
                self.quad = (Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3], Ssub(live[Q]), s3)

    def log_avoid(self, lam, order):
        P = T2.Pattern.__new__(T2.Pattern)
        P.SB = 1.0; P.s1 = self.s1; P.pairs = self.pairs; P.trip = self.trip; P.quad = self.quad
        if self.pairs is not None:
            P.s2 = self.s2
        return P.log_avoid(lam, order)


def theory():
    res = {}
    discs = {}
    for R in NORMS + [NORMS[-1] + 2]:
        t0 = time.time()
        d0 = Disc(R, 3)
        order = 4 if len(d0.s1) <= 40 else 3
        discs[R] = (Disc(R, order) if order == 4 else d0, order)
        print('disc R<%d live %d order %d %ds' % (R, len(d0.s1), order, time.time() - t0), flush=True)
    for box, (A, B0, K) in BOXES.items():
        res[box] = {}
        for blk in range(-1, 10):
            b0, b1 = (B0, B0 + K) if blk < 0 else (B0 + blk * K // 10, B0 + (blk + 1) * K // 10)
            Nmid = NRM(A + K / 2, (b0 + b1) / 2)
            lam_eff = math.log(Nmid) / KAPPA                  # p_w = rho_w * KAPPA / ln N
            ge = {R: math.exp(discs[R][0].log_avoid(lam_eff, discs[R][1])) for R in discs}   # P(d >= R)
            Rs = sorted(ge)
            pd = {Rs[i]: ge[Rs[i]] - ge[Rs[i + 1]] for i in range(len(Rs) - 1)}
            t = {('d=%d' % d): pd[d] for d in REPORT_D}
            t['d>=19'] = ge[19]
            t['mean_d'] = sum(R * v for R, v in pd.items()) + ge[Rs[-1]] * Rs[-1]
            res[box][str(blk)] = t
        print(box, {k: round(v, 5) for k, v in res[box]['-1'].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'theory_r16.json'), 'w'), indent=1)


def data():
    offs = sorted([(x, y) for x in range(-12, 13) for y in range(-12, 13) if 0 < NRM(x, y) <= RMAX],
                  key=lambda v: NRM(*v))
    rows = []
    for box, (A, B0, K) in BOXES.items():
        lo = A - MARGIN; loB = B0 - MARGIN; n = K + 2 * MARGIN
        M = np.zeros((n, n), dtype=bool)
        for i in range(n):
            a = lo + i
            M[i] = [gmpy2.is_prime(a * a - a * (loB + j) + (loB + j) ** 2) for j in range(n)]
        core = M[MARGIN:MARGIN + K, MARGIN:MARGIN + K]
        nn = np.full((K, K), 10 ** 6, dtype=np.int64)
        for (x, y) in offs:
            sh = M[MARGIN + x:MARGIN + x + K, MARGIN + y:MARGIN + y + K]
            upd = core & sh & (nn == 10 ** 6)
            nn[upd] = NRM(x, y)
        assert not np.any(core & (nn == 10 ** 6)), 'neighbour beyond RMAX'
        for blk in range(-1, 10):
            sl = slice(None) if blk < 0 else slice(blk * K // 10, (blk + 1) * K // 10)
            d = nn[:, sl][core[:, sl]]
            st = {'count': int(len(d))}
            for v in REPORT_D:
                st['d=%d' % v] = float(np.mean(d == v))
            st['d>=19'] = float(np.mean(d >= 19)); st['mean_d'] = float(np.mean(d))
            st['max_d'] = int(d.max())
            rows.append({'window': box, 'block': blk, 'stats': st})
        print(box, int(core.sum()), 'Eisenstein primes', flush=True)
    with open(os.path.join(OUT, 'results_r16.jsonl'), 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')


if __name__ == '__main__':
    {'data': data, 'theory': theory}[sys.argv[1]]()
