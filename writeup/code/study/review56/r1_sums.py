"""r1_sums.py (review, own code): prime sums for A4, B6, C23, C47a on the fresh chunks (3e7, XMAX],
chunk-matched against sums_fresh.jsonl; plus axis zone sums (m = 0, 4) for the Gaussian primes.
Independent methods: own segmented sieve over all integers; Q(i) and Q(sqrt-3) by enumeration over a (outer) with
lambda from the algebraic power (a+ib)^4 / (a+i sqrt3 b)^6; C23 by the Frobenius test x^p = x mod (x^3-x-1, p);
C47a by PARI qfbred(qfbprimeform(-47,p)).
usage: r1_sums.py XMAX"""
import os, sys, json, math, time
import numpy as np

SD = os.path.dirname(os.path.abspath(__file__))+'/'
FRESH = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/resonance/primeside/sums_fresh.jsonl'
try:
    XMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 193000000
    XMIN = int(float(sys.argv[2])) if len(sys.argv) > 2 else 0
    TAG = sys.argv[3] if len(sys.argv) > 3 else ''
except ValueError:
    XMAX, XMIN, TAG = 193000000, 0, ''
TP = 2 * math.pi
V0, V1 = math.sqrt(3e7), math.sqrt(1e9)          # registered main window (v = sqrt p)
CTRL = [0.6713, 0.7919, 0.8837, 0.9431, 1.0617, 1.1237, 1.2389, 1.3711]

def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)
BASE = small_primes(40000)

def prime_flags(L, U):
    """flags[k] True iff L+1+k is prime, k = 0 .. U-L-2 (numbers L+1 .. U-1)."""
    n = U - L - 1; s = np.ones(n, bool)
    for q in BASE:
        q = int(q)
        if q * q > U - 1: break
        st = max(q * q, -(-(L + 1) // q) * q)
        s[st - (L + 1)::q] = False
    return s

def isq(x):
    x = np.asarray(x, np.int64); r = np.sqrt(np.maximum(x, 0).astype(np.float64)).astype(np.int64)
    for _ in range(2):
        r = np.where(r * r > x, r - 1, r); r = np.where((r + 1) * (r + 1) <= x, r + 1, r)
    return r

def lattice_block(A, blo, bhi):
    """for each a in A: b in [blo, bhi] -> flat arrays."""
    cnt = np.maximum(bhi - blo + 1, 0); tot = int(cnt.sum())
    if tot == 0: return np.zeros(0, np.int64), np.zeros(0, np.int64)
    a = np.repeat(A, cnt); off = np.arange(tot, dtype=np.int64) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    return a, np.repeat(blo, cnt) + off

def gauss_primes(L, U, fl):
    """p = a^2 + b^2 prime in (L, U), a > b >= 1, enumerated over a (outer)."""
    amin = isq((L + 1) // 2); amax = isq(U - 1)
    outp, outa, outb = [], [], []
    for A0 in range(int(amin), int(amax) + 1, 2000):
        A = np.arange(A0, min(A0 + 2000, int(amax) + 1), dtype=np.int64)
        lo = np.where(L + 1 - A * A > 0, isq(np.maximum(L - A * A, 0)) + 1, 1)   # b^2 >= L+1-a^2
        lo = np.maximum(lo, 1); hi = np.minimum(A - 1, isq(U - 1 - A * A))
        a, b = lattice_block(A, lo, hi)
        n = a * a + b * b; m = (n > L) & (n < U)
        a, b, n = a[m], b[m], n[m]; m = fl[n - (L + 1)]
        outp.append(n[m]); outa.append(a[m]); outb.append(b[m])
    return np.concatenate(outp), np.concatenate(outa), np.concatenate(outb)

def eis_primes(L, U, fl):
    """p = a^2 + 3 b^2 prime, a >= 1, b >= 1, enumerated over a (outer)."""
    amax = isq(U - 1 - 3)
    outp, outa, outb = [], [], []
    for A0 in range(1, int(amax) + 1, 2000):
        A = np.arange(A0, min(A0 + 2000, int(amax) + 1), dtype=np.int64)
        r = L - A * A          # need 3b^2 > r
        lo = np.where(r >= 0, isq(np.maximum(r // 3, 0)) , 1); lo = np.maximum(lo, 1)
        hi = isq(np.maximum(U - 1 - A * A, 0) // 3)
        a, b = lattice_block(A, lo, hi)
        n = a * a + 3 * b * b; m = (n > L) & (n < U)
        a, b, n = a[m], b[m], n[m]; m = fl[n - (L + 1)]
        outp.append(n[m]); outa.append(a[m]); outb.append(b[m])
    return np.concatenate(outp), np.concatenate(outa), np.concatenate(outb)

def frob_split(p):
    """True iff x^3 - x - 1 splits completely mod p (x^p == x mod (f, p)); vectorised, p < 2^31."""
    p = p.astype(np.int64)
    def mul(A, B):
        a0, a1, a2 = A; b0, b1, b2 = B
        d0 = (a0 * b0) % p; d1 = ((a0 * b1) % p + (a1 * b0) % p) % p
        d2 = ((a0 * b2) % p + (a1 * b1) % p + (a2 * b0) % p) % p
        d3 = ((a1 * b2) % p + (a2 * b1) % p) % p; d4 = (a2 * b2) % p
        # x^3 = x + 1, x^4 = x^2 + x
        return ((d0 + d3) % p, (d1 + d3 + d4) % p, (d2 + d4) % p)
    def mulx(A):
        a0, a1, a2 = A; return (a2 % p, (a0 + a2) % p, a1 % p)
    R = (np.ones_like(p), np.zeros_like(p), np.zeros_like(p))
    nb = int(p.max()).bit_length()
    for k in range(nb - 1, -1, -1):
        R = mul(R, R); bit = ((p >> k) & 1).astype(bool)
        Rx = mulx(R); R = tuple(np.where(bit, x, y) for x, y in zip(Rx, R))
    return (R[0] == 0) & (R[1] == 1) & (R[2] == 0)

def leg(p, q):
    t = np.full(q, -1, np.int64); t[(np.arange(1, q) ** 2) % q] = 1; t[0] = 0
    return t[p % q]

_P = [None]
def pari():
    if _P[0] is None:
        import cypari2
        _P[0] = cypari2.Pari(); _P[0].allocatemem(200 * 10 ** 6)
    return _P[0]

C5a, C5b = 2 * math.cos(2 * math.pi / 5), 2 * math.cos(4 * math.pi / 5)
def c47_lam(p):
    """lambda for C47a: 2 (principal), 2cos(2pi/5) (class of (2,1,6)), 2cos(4pi/5) (class of (3,1,4)), 0 inert."""
    lam = np.zeros(len(p)); sp = leg(p, 47) == 1
    ps = p[sp]
    P = pari(); out = []
    for i in range(0, len(ps), 200000):
        v = P('my(v=%s); Vecsmall(vector(#v, i, component(qfbred(qfbprimeform(-47, v[i])), 1)))' % ('[' + ','.join(map(str, ps[i:i + 200000].tolist())) + ']'))
        out.append(np.array(list(v), np.int64))
    a = np.concatenate(out) if out else np.zeros(0, np.int64)
    assert set(np.unique(a).tolist()) <= {1, 2, 3}
    lam[sp] = np.where(a == 1, 2.0, np.where(a == 2, C5a, C5b))
    return lam

def sums(p, lam, Q, Vw=None):
    pf = p.astype(np.float64); lp = np.log(pf); v = np.sqrt(pf)
    ph = TP * 2 * np.sqrt(pf / Q)
    u = (v - V0) / (V1 - V0); W = np.where((u > 0) & (u < 1), np.sin(np.pi * u) ** 2, 0.0)
    g = pf ** -0.25; wl = lam * lp; ws = wl * g * W
    E = np.exp(-1j * ph)
    out = {'sharp': complex(np.dot(wl, E)), 'smooth': complex(np.dot(ws, E)),
           'noise_sm': float(np.dot(ws, ws)), 'ctrl_sm': [complex(np.dot(ws, np.exp(-1j * a * ph))) for a in CTRL]}
    if Vw is not None:   # own window on the sub-range
        u2 = (v - Vw[0]) / (Vw[1] - Vw[0]); W2 = np.where((u2 > 0) & (u2 < 1), np.sin(np.pi * u2) ** 2, 0.0)
        w2 = wl * g * W2
        out['mine'] = complex(np.dot(w2, E)); out['noise_mine'] = float(np.dot(w2, w2))
        out['ctrl_mine'] = [complex(np.dot(w2, np.exp(-1j * a * ph))) for a in CTRL]
    return out

def zones(p, a, b):
    pf = p.astype(np.float64); af = a.astype(np.float64); bf = b.astype(np.float64)
    lp = np.log(pf); base = lp * np.exp(-1j * TP * np.sqrt(pf))
    z = (af + 1j * bf) / np.sqrt(pf); e4 = z ** 4
    u = bf / np.sqrt(af); w = np.exp(-np.pi * u * u / 4); w[u > 12] = 0
    return {'ax0_m0': complex(np.sum(w * base)), 'ax0_m4': complex(np.sum(w * base * e4)),
            'ax90_m4': complex(np.sum(w * base * np.conj(e4))), 'noise': float(np.sum((w * lp) ** 2)),
            # the same with a sharp angular cut (no window): diagnostic
            }

def main():
    recs = [json.loads(l) for l in open(FRESH)]
    abc = sorted([(r['L'], r['U'], r) for r in recs if r['group'] == 'ABC'], key=lambda t: t[0])
    abc = [t for t in abc if t[1] <= XMAX and t[0] >= XMIN]
    U_end = abc[-1][1]; Vw = (V0, math.sqrt(193263832))
    print('chunks', len(abc), 'from', abc[0][0], 'to', U_end, flush=True)
    res = []; t0 = time.time()
    for (L, U, r) in abc:
        t1 = time.time(); fl = prime_flags(L, U)
        allp = L + 1 + np.flatnonzero(fl)
        out = {'L': L, 'U': U, 'nprimes': int(len(allp))}
        pA, aA, bA = gauss_primes(L, U, fl)
        assert np.all(pA % 4 == 1) and len(np.unique(pA)) == len(pA) == int(np.sum(allp % 4 == 1)), 'Q(i) enumeration'
        za = (aA.astype(np.float64) + 1j * bA) ; lamA = 2 * np.real((za / np.abs(za)) ** 4)
        out['A4'] = sums(pA, lamA, 4, Vw); out['zones'] = zones(pA, aA, bA)
        pB, aB, bB = eis_primes(L, U, fl)
        assert np.all(pB % 3 == 1) and len(np.unique(pB)) == len(pB) == int(np.sum(allp % 3 == 1)), 'Q(sqrt-3) enumeration'
        zb = (aB.astype(np.float64) + 1j * math.sqrt(3) * bB); lamB = 2 * np.real((zb / np.abs(zb)) ** 6)
        out['B6'] = sums(pB, lamB, 3, Vw)
        lg = leg(allp, 23); spl = frob_split(allp)
        assert not np.any(spl & (lg != 1)), 'split but (p/23) != 1'
        lamC = np.where(lg == 1, np.where(spl, 2.0, -1.0), 0.0)
        out['C23'] = sums(allp, lamC, 23, Vw)
        lam47 = c47_lam(allp)
        out['C47a'] = sums(allp, lam47, 47, Vw)
        out['counts'] = {'A4': int(len(pA)), 'B6': int(len(pB)), 'C23_split': int(spl.sum()), 'C47_principal': int(np.sum(lam47 == 2.0))}
        # comparison with the file
        cmp = {}
        for k in ('A4', 'B6', 'C23', 'C47a'):
            q = r['res'][k]; fs = complex(*q['sharp'][0]); fw = complex(*q['smooth'][0][0])
            cmp[k] = [abs(out[k]['sharp'] - fs) / max(1.0, abs(fs)), abs(out[k]['smooth'] - fw) / max(1.0, abs(fw))]
        for nm, key in (('axis0', 'ax0_m0'), ('axis0', 'ax0_m4'), ('axis90', 'ax90_m4')):
            fz = complex(*r['zones'][nm][key[-1]]); cmp[nm + key[-3:]] = abs(out['zones'][key] - fz) / max(1.0, abs(fz))
        out['cmp'] = cmp
        res.append(out)
        print('(%d,%d] %d primes  %.1fs  max rel diff %.2e' % (L, U, len(allp), time.time() - t1, max(max(v) if isinstance(v, list) else v for v in cmp.values())), flush=True)
    def enc(o):
        if isinstance(o, complex): return [o.real, o.imag]
        raise TypeError(type(o))
    json.dump({'Vw': Vw, 'chunks': res}, open(SD + '/r1_out%s.json' % TAG, 'w'), default=enc)
    print('total %.0fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
