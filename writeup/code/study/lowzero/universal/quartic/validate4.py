# independent checks: (1) per-field bin sums recomputed with plain integers (x^p mod f mod p, gcd), sympy primes; (2) the Mellin evaluation
# reproduces the S3 cubic-field table of round 161 away from its spikes; (3) Dirichlet coefficients B(rho) from the closed Euler factor
# against the direct type average; small-y coefficients.
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sympy as sp
H = os.path.dirname(os.path.abspath(__file__))
def pmod(a, f, p):          # a mod (f, p), f monic, lists low->high
    a = [c % p for c in a]; d = len(f) - 1
    while len(a) > d:
        c = a.pop()
        if c:
            for i in range(d): a[len(a) - d + i] = (a[len(a) - d + i] - c * f[i]) % p
    return a
def pmul(a, b, f, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] += x * y
    return pmod(r, f, p)
def pgcd_deg(a, b, p):
    def trim(v):
        while v and v[-1] % p == 0: v.pop()
        return v
    a, b = trim([c % p for c in a]), trim([c % p for c in b])
    while b:
        inv = pow(b[-1], -1, p)
        while len(a) >= len(b) and a:
            c = a[-1] * inv % p; s = len(a) - len(b)
            for i in range(len(b)): a[s + i] = (a[s + i] - c * b[i]) % p
            a = trim(a)
        a, b = b, a
    return len(a) - 1
def nroots(f, p):
    r = [1]; base = [0, 1]; e = p
    while e:
        if e & 1: r = pmul(r, base, f, p)
        base = pmul(base, base, f, p); e >>= 1
    r = r + [0] * (4 - len(r)); r[1] = (r[1] - 1) % p        # x^p - x mod f
    return pgcd_deg(f[:], r, p)
x = sp.symbols('x'); FE = np.linspace(0, 2.0, 193); worst = 0.0; nchk = 0
for r2 in (1, 2):
    z = np.load(os.path.join(H, 'data', 'disc_r%d.npz' % r2)); rows = [l.split() for l in open(os.path.join(H, 'data', 'disc_r%d_fields.txt' % r2)) if not l.startswith('#')]
    rows = sorted(set((int(D), f) for D, f in rows), key=lambda r: (abs(r[0]), r[1])); picked = 0
    for i, (D, fs) in enumerate(rows):
        poly = sp.Poly(sp.sympify(fs.replace('^', '**')), x)
        if picked >= 3: continue
        ind = math.isqrt(int(sp.discriminant(poly)) // D)
        if ind * ind * D != int(sp.discriminant(poly)) or any(D % q for q in sp.primefactors(ind)): continue   # index primes must divide D (then dropped)
        f = [int(c) for c in reversed(poly.all_coeffs())]; A = abs(D); num = np.zeros(192); den = np.zeros(192)
        for p in sp.primerange(2, 8 * A + 1):
            if A % p == 0: continue
            u = (p / A) ** (1 / 3); j = int(np.searchsorted(FE, u, side='left')) - 1
            if 0 <= j < 192: a = nroots(f, p) - 1; num[j] += math.log(p) * math.sqrt(A) * a; den[j] += math.log(p)
        d = max(np.abs(num - z['num'][i]).max(), np.abs(den - z['den'][i]).max()); worst = max(worst, d); picked += 1; nchk += 1
        print('r2=%d field D=%d %s: max |difference| of the 192 bin sums %.2e' % (r2, D, fs, d), flush=True)
print('(1) %d fields recomputed from scratch, largest difference %.1e' % (nchk, worst))
if '--fields-only' in sys.argv: sys.exit()
import pred4 as P
tab = np.loadtxt('/home/bas/praime/writeup/code/study/lowzero/universal/cubic/m_table.txt'); cub = P.cubic_model(T=2000.0)
spikes = [k * k / 4 for k in range(1, 4)] + [1 / 8, 9 / 8, 1 / 12, 1 / 3, 4 / 3, 25 / 12, 1 / 16, 9 / 16] + [(2 * k + 1)**2 / 32 for k in range(5)]
e = np.round(np.arange(0.01, 2.99, 0.03), 10)
keep = np.array([min(abs(s - a) for s in spikes) > 0.012 and min(abs(s - b) for s in spikes) > 0.012 for a, b in zip(e[:-1], e[1:])])
for sg, col in ((-1, 2), (1, 3)):
    ref = np.array([tab[int(round(a * 100)):int(round(b * 100)), col].mean() for a, b in zip(e[:-1], e[1:])]); d = P.bin_means(cub, P.CUBIC_MU[sg], e) - ref
    print('(2) cubic D%s0: %d of %d bins of width 0.03 with no edge within 0.012 of a listed spike: max |Mellin - table| %.4f, rms %.4f' % ('<' if sg < 0 else '>', keep.sum(), len(keep), np.abs(d[keep]).max(), np.sqrt(np.mean(d[keep]**2))))
# (3) B(rho) for rho = l^v from the type average and from the closed Euler factor, and the y^(-1/4) coefficient
for l in (2, 3, 7):
    Xs = np.array([0.11, 0.23]); a = P.E_local(l, Xs); b = (1 + l**-2 - l**-3 - l**-4)**-1 * P.P_l(l, Xs) / ((1 - Xs**2) * (1 - Xs**3) * (1 - Xs**4))
    print('(3) l=%d: E_l(X) by types %s, closed form %s' % (l, np.round(a, 12), np.round(b, 12)))
A1 = P.A1_exact(); q4 = P.Q_real_acc(0.25); r4 = A1 * P.zeta_real(0.5) * P.zeta_real(0.75) * q4 / 4
for r2, mus in P.SIG_MU.items(): print('    r2=%d: y^(-1/4) coefficient of m: %.4f' % (r2, r4 * float(P.Kt(np.array([0.25]), mus)[0].real)))
