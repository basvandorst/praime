"""Independent check of the high-precision chain (hpchain.py): mpmath and sympy, no Arb.

For each recorded window the script rebuilds everything itself:
- the zeros with mpmath.zetazero (not Arb's), at 50 digits;
- each sample: zero side with mpmath's Bessel functions, trivial zeros by mp.quad, prime side exact over the
  numbers outside the window (true Lambda from sympy), columns for the unknowns;
- the implied rows from the recorded weights u and scale factors: sum_r u_r (A b)_r lies between
  sum u+ lo - u- hi and sum u+ hi - u- lo;
- delta: only the zeros above T (cert.tail_bound, cert.offline_bound), as in hpwin.py.
It then replays every recorded proof tree: forcings and violated rows on the sample rows, Farkas vectors on all rows,
and the tree structure, with comparisons at 50 digits. It also checks the identity at random sample points.
Usage: hp_check.py run.json [max_windows]   (the trees are read from run_proofs.json or run_proofs.json.gz)
"""
import json, math, random, sys, os, gzip
import mpmath as mp
from sympy import isprime, factorint, primerange
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cert import tail_bound, offline_bound

mp.mp.dps = 50
_Z = []
def zeros_upto(T):
    while not _Z or _Z[-1] <= T:
        _Z.append(mp.zetazero(len(_Z) + 1).imag)
    return [g for g in _Z if g <= T]

# a, beta and x arrive as Python floats: convert them before any arithmetic, or beta ** 2 and a / (2 x) are rounded
# to double precision (an error of about 1e-13 K0 at beta = 40, larger than delta)
def Kmp(v, a, beta, m):
    a = mp.mpf(a); beta = mp.mpf(beta)
    t = 1 - (v / a) ** 2
    return t ** (mp.mpf(m) / 2) * mp.besseli(m, beta * mp.sqrt(t)) if t > 0 else mp.mpf(0)

def Khmp(w, a, beta, m):
    a = mp.mpf(a); beta = mp.mpf(beta)
    nu = m + mp.mpf(1) / 2; s = abs(a * w)
    if s < beta:
        q = mp.sqrt(beta ** 2 - s ** 2); z = mp.besseli(nu, q) / q ** nu if q > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    else:
        y = mp.sqrt(s ** 2 - beta ** 2); z = mp.besselj(nu, y) / y ** nu if y > 0 else 1 / (2 ** nu * mp.gamma(nu + 1))
    return a * mp.sqrt(2 * mp.pi) * beta ** m * z

def Pmp(x, a, beta, m):
    x = mp.mpf(x); a = mp.mpf(a); beta = mp.mpf(beta)
    nu = m + mp.mpf(1) / 2; r = a / (2 * x); q = mp.sqrt(beta ** 2 + r ** 2)
    return a * mp.sqrt(2 * mp.pi) * beta ** m * mp.besseli(nu, q) / q ** nu

_lam = {}
def lam(n):
    if n not in _lam:
        f = factorint(n); _lam[n] = mp.log(next(iter(f))) if len(f) == 1 else mp.mpf(0)
    return _lam[n]

def zero_side(x, a, beta, m, T):
    x = mp.mpf(x); L = mp.log(x); s = mp.mpf(0)
    for g in zeros_upto(T):
        s += mp.cos(g * L) * Khmp(g / x, a, beta, m)
    triv = mp.quad(lambda v: Kmp(v, a, beta, m) * mp.sqrt(x * mp.e ** (v / x)) / (x * ((x * mp.e ** (v / x)) ** 3 - x * mp.e ** (v / x))), [-a, 0, a])
    return Pmp(x, a, beta, m) - 2 / mp.sqrt(x) * s - mp.sqrt(x) * triv

def delta(x, a, beta, m, T):
    NT = len(zeros_upto(T))
    return ((2 / math.sqrt(x)) * tail_bound(x, T, NT, a, beta, m) + offline_bound(x, a, beta, m)) * (1 + 1e-6)

def support(x, a):
    x = mp.mpf(x); return int(mp.floor(x * mp.e ** (-a / x))) - 1, int(mp.ceil(x * mp.e ** (a / x))) + 1

out = []
def say(s):
    print(s, flush=True); out.append(s)

J = json.load(open(sys.argv[1])); maxwin = int(sys.argv[2]) if len(sys.argv) > 2 else 1
start = J.get('start', 3); X = J['X']; hi_dec = J['result']['stop'] or X
true = [p for p in primerange(max(3, start), hi_dec)]
got = [p for p in J['primes'] if p < hi_dec]
say(f'1. {len(got)} certified primes in [{max(3, start)}, {hi_dec}); equal to sympy: {got == true}; stop: {J["result"]["stop"]}; replays in the run: {J["result"].get("replays")}')

random.seed(11)
for st in J['stages'][-2:]:
    end, a, beta, m, cc, k = st
    x = float(random.randint(max(int(a * 1.5), end // 2), end - 1)); T = cc * beta * x / a
    lo_, hi_ = support(x, a); xm = mp.mpf(x)
    ps = sum(lam(n) * mp.sqrt(xm / n) * Kmp(xm * mp.log(mp.mpf(n) / xm), a, beta, m) for n in range(max(2, lo_), hi_ + 1))
    zs = zero_side(x, a, beta, m, T); d = delta(x, a, beta, m, T); k0 = mp.besseli(m, beta)
    say(f'2. x = {x:.0f} (a = {a}, beta = {beta}, T = {T:.1f}): |zero side - prime side| / K0 = {mp.nstr(abs(zs - ps) / k0, 3)}, delta / K0 = {d / float(k0):.3e}, inside: {abs(zs - ps) <= d}')

pf = sys.argv[1].replace('.json', '_proofs.json')
P = json.load(open(pf)) if os.path.exists(pf) else json.load(gzip.open(pf + '.gz', 'rt')) if os.path.exists(pf + '.gz') else {}
if not P: say(f'3. no recorded proofs ({os.path.basename(pf)}): window replay skipped')
keys = sorted(P, key=lambda c: -sum(len(p['nodes']) for p in P[c]['proofs']))[:maxwin]
for key in keys:
    R = P[key]; U = R['U']; xs = R['xs']; a, beta, m, cc, k, h = R['params']; n0 = R['n0']
    k0 = mp.besseli(m, beta); Uset = set(U)
    A = []; lo = []; hi = []
    for x in xs:
        T = cc * beta * x / a; xm = mp.mpf(x)
        l_, h_ = support(x, a)
        known = sum(lam(n) * mp.sqrt(xm / n) * Kmp(xm * mp.log(mp.mpf(n) / xm), a, beta, m) for n in range(max(2, l_), h_ + 1) if n not in Uset)
        y = (zero_side(x, a, beta, m, T) - known) / k0; d = mp.mpf(delta(x, a, beta, m, T)) / k0
        lo.append(y - d); hi.append(y + d)
        A.append([mp.log(u) * mp.sqrt(xm / u) * Kmp(xm * mp.log(mp.mpf(u) / xm), a, beta, m) / k0 for u in U])
    for wv, sc in zip(R['weights'], R['scales']):
        row = [sum(mp.mpf(wv[r]) * A[r][j] for r in range(n0)) * sc for j in range(len(U))]
        l = sum((mp.mpf(wv[r]) * (lo[r] if wv[r] > 0 else hi[r]) for r in range(n0))) * sc
        hh = sum((mp.mpf(wv[r]) * (hi[r] if wv[r] > 0 else lo[r]) for r in range(n0))) * sc
        A.append(row); lo.append(l); hi.append(hh)
    nr = len(A); kk = len(U)
    def rows(i, st):
        base = mp.mpf(0); tot = mp.mpf(0)
        for j in range(kk):
            if st[j] == 1: base += A[i][j]
            elif st[j] == -1: tot += A[i][j]
        return base, tot
    bad = []; certified = {}
    for pr in R['proofs']:
        j, val = pr['bit'], pr['value']
        dec = lambda s: [{'x': -1, '0': 0, '1': 1}[ch] for ch in s]
        r0 = dec(pr['root'])
        if r0[j] != 1 - val: bad.append(f'bit {j}: root does not set 1 - value')
        if any(r0[i] != certified[i] for i in certified if i != j) or any(r0[i] != -1 for i in range(kk) if i != j and i not in certified):
            bad.append(f'bit {j}: root state does not fix exactly the bits certified before it')
        states = set(); need = []
        for node in pr['nodes']:
            st = dec(node['state']); states.add(node['state']); closed = False
            for ev in node['events']:
                if ev[0] == 'force':
                    _, var, v, row = ev
                    base, tot = rows(row, st)
                    ok = row < n0 and (base + A[row][var] > hi[row] if v == 0 else base + tot - A[row][var] < lo[row])
                    if st[var] != -1 or not ok: bad.append(f'bit {j}: forcing of {U[var]} by row {row} fails')
                    st[var] = v
                else:
                    _, row, side = ev
                    base, tot = rows(row, st)
                    ok = row < n0 and (base > hi[row] if side == 'hi' else base + tot < lo[row])
                    if not ok: bad.append(f'bit {j}: row {row} is not violated')
                    closed = True
            if closed: continue
            if node.get('farkas'):
                w = [mp.mpf(v) for v in node['farkas']]; free = [q for q in range(kk) if st[q] == -1]
                nz = [i for i in range(nr) if w[i] != 0]
                g = [sum(w[i] * A[i][q] for i in nz) for q in free]
                minbox = sum(min(gq, 0) for gq in g)
                rhs = sum(w[i] * ((hi[i] if w[i] > 0 else lo[i]) - rows(i, st)[0]) for i in nz)
                if not minbox > rhs: bad.append(f'bit {j}: Farkas vector fails ({mp.nstr(minbox - rhs, 3)})')
            elif 'branch' in node:
                q = node['branch']
                for v in (0, 1):
                    s2 = list(st); s2[q] = v; need.append(''.join('x01'[t + 1] for t in s2))
            else:
                bad.append(f'bit {j}: open node')
        if any(s not in states for s in need): bad.append(f'bit {j}: children missing')
        certified[j] = val
    agree = all(certified[j] == (1 if isprime(U[j]) else 0) for j in certified)
    say(f'3. window at {U[0]} (a = {a}, beta = {beta}, {n0} samples + {nr - n0} implied rows, {kk} unknowns): {len(certified)} decisions replayed, '
        f'{sum(len(p["nodes"]) for p in R["proofs"])} nodes, failures: {bad[:3] if bad else 0}; decisions equal sympy: {agree}')
open(sys.argv[1].replace('.json', '_check.txt'), 'w').write('\n'.join(out) + '\n')
