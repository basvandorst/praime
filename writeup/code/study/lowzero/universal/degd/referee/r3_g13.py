"""Referee check 3: G1 (zeros = primes + E0, Weil form) and G3 / G3' (dual side), own data and code.
Cases: chi mod 7 (order 6, odd; complex eps), Sym^2 11a (degree 3).
Zeros: PARI lfunzeros. Coefficients: chi from PARI chareval; Sym^2 from ellan(E) (own Satake algebra, own local factor at 11),
cross-checked against PARI lfunan(lfunsympow). Lambda_L(p^k) = (sum_j alpha_j^k) ln p via Chebyshev polynomials.
Quadrature: trapezoid + Horner (different from the note's restart recurrence)."""
import sys, math, cmath, time
import numpy as np
from scipy.special import loggamma, digamma, erf
import cypari2
P = cypari2.Pari(); P.allocatemem(800 * 10**6); P.set_real_precision(38)
PI = math.pi

def lnGR(z): return -(z / 2) * math.log(PI) + loggamma(z / 2)
def e2it(t, Q, mu):   # X^{-1}(1/2+it) directly
    s = 0.5 + 1j * np.asarray(t, float)
    v = 1j * np.asarray(t, float) * math.log(Q) + 0j
    for m in mu: v = v + lnGR(s + m) - lnGR(1 - s + np.conj(m))
    return np.exp(v)
def thetap(t, Q, mu):  # theta'(t) = (1/2)(-X'/X)(1/2+it)
    s = 0.5 + 1j * np.asarray(t, float)
    v = math.log(Q) + 0j
    for m in mu: v = v + (-math.log(PI) + 0.5 * digamma((s + m) / 2)) * 0.5 * 2 / 2 * 2 / 2 * 2 + 0  # placeholder replaced below
    return None
def psiR(z): return -0.5 * math.log(PI) + 0.5 * digamma(z / 2)
def thetap(t, Q, mu):
    s = 0.5 + 1j * np.asarray(t, float)
    v = math.log(Q) + 0j
    for m in mu: v = v + psiR(s + m) + psiR(1 - s + np.conj(m))
    return 0.5 * v.real, np.abs(v.imag).max()
def Hbox(t, T1, T2, D): return 0.5 * (erf((np.asarray(t) - T1) / D) - erf((np.asarray(t) - T2) / D))
def Hhat(xi, T1, T2, D):
    xi = np.asarray(xi, float)
    W = np.where(np.abs(xi) < 1e-14, T2 - T1, (np.exp(-1j * xi * T1) - np.exp(-1j * xi * T2)) / (1j * np.where(np.abs(xi) < 1e-14, 1, xi)))
    return W * np.exp(-D * D * xi * xi / 4)
def horner(L, t0, h, f, sign=-1):
    """h * sum_k f_k exp(sign i (t0 + k h) L), vectorised over L."""
    z = np.exp(sign * 1j * h * L)
    acc = np.full(L.shape, f[-1], complex)
    for k in range(len(f) - 2, -1, -1):
        acc = acc * z + f[k]
    return h * np.exp(sign * 1j * t0 * L) * acc
def primes(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

out = []
def say(x): print(x, flush=True); out.append(x)

# ---------------------------------------------------------------- chi mod 7
G = P('znstar(7,1)'); L7 = P.lfuncreate(P('[znstar(7,1),[1]]')); L7b = P.lfuncreate(P('[znstar(7,1),[5]]'))
eps7 = complex(P.lfunrootres(L7)[2]); Q7, mu7 = 7.0, [1.0]
chi = np.array([0] + [complex(P(f'chareval(znstar(7,1),[1],{n},[exp(2*Pi*I/6),6])'.replace('[exp(2*Pi*I/6),6]', '[-1,0]')) if False else 0) for n in range(1, 8)])
# chareval with z = exp(2 pi i/6): value e(chi-log)
chi = np.zeros(7, complex)
for n in range(1, 7):
    v = P(f'chareval(znstar(7,1),[1],{n})')   # rational x, chi(n) = e(x)
    chi[n] = cmath.exp(2j * PI * float(v))
say(f'chi mod 7: values {np.round(chi[1:], 6)}, eps = {eps7:.12f}')
t0 = time.time()
z7 = np.array([float(x) for x in P.lfunzeros(L7, P('[0,290]'))])
z7b = np.array([float(x) for x in P.lfunzeros(L7b, P('[0,290]'))])
say(f'zeros chi7: {len(z7)} in [0,290]; chi7bar: {len(z7b)}; {time.time()-t0:.1f}s')
N = 40000
pr = primes(N)
ns, Lam = [], []
for p in pr:
    p = int(p); q = p; k = 1
    while q <= N:
        ns.append(q); Lam.append(chi[p % 7] ** k * math.log(p)); q *= p; k += 1
ns = np.array(ns, float); Lam = np.array(Lam)
for (T1, T2, D, r) in ((3.0, 15.0, 1.5, 1.0), (3.0, 15.0, 1.5, 2.5), (200.0, 230.0, 4.0, 1.0), (200.0, 230.0, 4.0, 2.5), (190.0, 240.0, 6.0, 1.0)):
    h = 0.05; a, b = T1 - 6.5 * D, T2 + 6.5 * D
    t = a + h * np.arange(int((b - a) / h) + 1)
    f = Hbox(t, T1, T2, D) * e2it(t, Q7, mu7)
    gpos = z7[(z7 > a - 1) & (z7 < b + 1)]; gneg = -z7b[(-z7b > a - 1) & (-z7b < b + 1)]  # zeros of L at -gamma(Lbar)
    g = np.concatenate([gpos, gneg])
    Z = np.sum(Hbox(g, T1, T2, D) * e2it(g, Q7, mu7) * np.exp(-1j * g * math.log(r)))
    I = horner(np.log(ns) + math.log(r), t[0], h, f)
    Pside = -np.sum(Lam * ns ** -0.5 * I) / (2 * PI)
    Pconj = -np.sum(np.conj(Lam) * ns ** -0.5 * I) / (2 * PI)
    tail = abs(np.sum((Lam * ns ** -0.5 * I)[ns > N / 2])) / (2 * PI)
    thp, imres = thetap(t, Q7, mu7)
    Im_ = horner(np.log(ns) - math.log(r), t[0], h, f, sign=+1)
    E0 = h * np.sum(f * thp * np.exp(-1j * t * math.log(r))) / PI - np.sum(np.conj(Lam) * ns ** -0.5 * Im_) / (2 * PI)
    E0wrong = h * np.sum(f * thp * np.exp(-1j * t * math.log(r))) / PI - np.sum(Lam * ns ** -0.5 * Im_) / (2 * PI)
    say(f'G1 chi7 [{T1},{T2}] D={D} r={r}: zeros {Z:.10f} ({len(g)} zeros), primes {Pside:.10f}, |E0| {abs(E0):.1e}, '
        f'|zeros-primes-E0| = {abs(Z - Pside - E0):.2e}; tail {tail:.1e}; with conj Lambda on the resonant side: {abs(Z - Pconj - E0):.3f}; '
        f'E0 with unconjugated Lambda: {abs(Z - Pside - E0wrong):.2e}')
    # G3 (d = 1) for the same window, all n
    nn = np.arange(1, 6001, dtype=float)
    lam = chi[(nn.astype(int)) % 7]
    In = horner(np.log(nn) + math.log(r), t[0], h, f)
    lhs = np.sum(lam * nn ** -0.5 * In)
    m = np.arange(1, 4001, dtype=float); lm = chi[m.astype(int) % 7]
    dual = np.sum(np.conj(lm) * m ** -0.5 * Hhat(np.log(r / m), T1, T2, D))
    say(f'G3 chi7 r={r}: lhs {lhs:.10f}, eps*dual {eps7*dual:.10f}, |diff| {abs(lhs - eps7*dual):.2e}; '
        f'conj(eps): {abs(lhs - eps7.conjugate()*dual):.3f}; no conj on lambda: {abs(lhs - eps7*np.sum(lm * m ** -0.5 * Hhat(np.log(r / m), T1, T2, D))):.3f}')

# ---------------------------------------------------------------- Sym^2 11a
E = P.ellinit([0, -1, 1, -10, -20]); LS = P.lfunsympow(E, 2)
epsS = complex(P.lfunrootres(LS)[2]); QS, muS = 121.0, [1.0, 1.0, 2.0]
N = 700000
t0 = time.time()
an = np.array([int(x) for x in P.ellan(E, N)], dtype=float)  # a_1..a_N
say(f'ellan to {N}: {time.time()-t0:.1f}s; eps(Sym2 11a) = {epsS}')
pr = primes(N)
lam = np.ones(N + 1, complex); lam[0] = 0
ppn, ppL = [], []
for p in pr:
    p = int(p)
    K = int(math.log(N) / math.log(p)) + 1
    if p == 11:
        loc = [11.0 ** -k for k in range(K + 1)]
        sk = [11.0 ** -k for k in range(K + 1)]
    else:
        x = an[p - 1] / math.sqrt(p)           # beta + 1/beta
        A = x * x - 1                             # lambda_Sym2(p)
        # local factor 1/((1 - b^2 X)(1 - X)(1 - b^-2 X)) = 1/(1 - A X + A X^2 - X^3)
        loc = [1.0, A]
        for k in range(2, K + 1):
            loc.append(A * loc[k - 1] - A * loc[k - 2] + (loc[k - 3] if k >= 3 else 0.0))
        # power sums b^{2k} + 1 + b^{-2k} = 1 + 2 T_{2k}(x/2)
        Tm = [1.0, x / 2]
        for j in range(2, 2 * K + 1): Tm.append(2 * (x / 2) * Tm[-1] - Tm[-2])
        sk = [None] + [1 + 2 * Tm[2 * k] for k in range(1, K + 1)]
    q = p; k = 1
    while q <= N:
        idx = np.arange(q, N + 1, q)
        if q * p <= N: idx = idx[(idx // q) % p != 0]
        lam[idx] *= loc[k]
        ppn.append(q); ppL.append(sk[k] * math.log(p))
        q *= p; k += 1
ppn = np.array(ppn, float); ppL = np.array(ppL, complex); o = np.argsort(ppn); ppn, ppL = ppn[o], ppL[o]
anS = P.lfunan(LS, 3000)
ref = np.array([complex(anS[i]) for i in range(3000)]) / np.arange(1, 3001)
say(f'Sym2 11a coefficients (own, from ellan) vs PARI lfunan(lfunsympow), n<=3000: max diff {np.max(np.abs(ref - lam[1:3001])):.1e}; built in {time.time()-t0:.1f}s')
# Lambda check: k lam(p^k) ln p = sum_j Lambda(p^j) lam(p^{k-j}) on a few p
errL = 0
D_ = dict(zip(ppn.astype(np.int64).tolist(), ppL))
for p in (2, 3, 5, 7, 11, 13):
    for k in range(1, 6):
        lhs = k * lam[p ** k] * math.log(p); rhs = sum(D_[p ** j] * lam[p ** (k - j)] for j in range(1, k + 1))
        errL = max(errL, abs(lhs - rhs))
say(f'Lambda_L via Chebyshev vs Newton identity: max err {errL:.1e}')

t0 = time.time()
zS = np.array([float(x) for x in P.lfunzeros(LS, 75)])
say(f'zeros Sym2 11a up to 75: {len(zS)}, first {zS[:3]}; {time.time()-t0:.1f}s')
for (T1, T2, D, r) in ((30.0, 40.0, 5.0, 1.0), (30.0, 40.0, 5.0, 5 / 3)):
    h = 0.05; a, b = T1 - 6.5 * D, T2 + 6.5 * D
    t = a + h * np.arange(int((b - a) / h) + 1)
    f = Hbox(t, T1, T2, D) * e2it(t, QS, muS)
    g = np.concatenate([zS, -zS])
    Z = np.sum(Hbox(g, T1, T2, D) * e2it(g, QS, muS) * np.exp(-1j * g * math.log(r)))
    t1 = time.time()
    I = horner(np.log(ppn) + math.log(r), t[0], h, f)
    Pside = -np.sum(ppL * ppn ** -0.5 * I) / (2 * PI)
    tail = abs(np.sum((ppL * ppn ** -0.5 * I)[ppn > N / 2])) / (2 * PI)
    thp, imres = thetap(t, QS, muS)
    sel = ppn <= 5000
    Im_ = horner(np.log(ppn[sel]) - math.log(r), t[0], h, f, sign=+1)
    E0 = h * np.sum(f * thp * np.exp(-1j * t * math.log(r))) / PI - np.sum(np.conj(ppL[sel]) * ppn[sel] ** -0.5 * Im_) / (2 * PI)
    say(f'G1 Sym2 11a [{T1},{T2}] D={D} r={r:.4f}: zeros {Z:.10f} ({np.sum(np.abs(g) < 75)} zeros +-), primes {Pside:.10f}, |E0| {abs(E0):.1e}, '
        f'|zeros-primes-E0| = {abs(Z - Pside - E0):.2e}; tail(n>N/2) {tail:.1e}; zero-window margin H(75) = {Hbox(75.0, T1, T2, D):.1e}; {time.time()-t1:.1f}s')
# G3, G3' for Sym2 11a, all n <= N
T1, T2, D = 30.0, 40.0, 5.0
h = 0.05; a, b = T1 - 6.5 * D, T2 + 6.5 * D
t = a + h * np.arange(int((b - a) / h) + 1)
Ht = Hbox(t, T1, T2, D); f = Ht * e2it(t, QS, muS)
nn = np.arange(1, N + 1, dtype=float)
thp, _ = thetap(t, QS, muS)
for r in (1.5, 5 / 3, 3.0):
    t1 = time.time()
    In = horner(np.log(nn) + math.log(r), t[0], h, f)
    terms = lam[1:] * nn ** -0.5 * In
    lhs = np.sum(terms); tail = abs(np.sum(terms[nn > N / 2]))
    m = np.arange(1, 3001, dtype=float)
    dual = epsS * np.sum(np.conj(lam[1:3001]) * m ** -0.5 * Hhat(np.log(r / m), T1, T2, D))
    say(f'G3 Sym2 11a r={r:.4f}: lhs {lhs:.10f}, eps*dual {dual:.10f}, |diff| {abs(lhs - dual):.2e}, flipped eps {abs(lhs + dual):.3f}, tail {tail:.1e}; {time.time()-t1:.1f}s')
    # G3': sum lam(b) ln b b^{-1/2} I_b(a) = eps sum conj lam(m) m^{-1/2} int H (2 theta' - ln m) e^{-it ln(a/m)}
    lhs2 = np.sum(lam[1:] * np.log(nn) * nn ** -0.5 * In)
    rhs2 = 0j
    for mm in range(1, 1500):
        rhs2 += np.conj(lam[mm]) * mm ** -0.5 * h * np.sum(Ht * (2 * thp - math.log(mm)) * np.exp(-1j * t * math.log(r / mm)))
    rhs2 *= epsS
    say(f"G3' Sym2 11a a={r:.4f}: lhs {lhs2:.10f}, dual {rhs2:.10f}, |diff| {abs(lhs2 - rhs2):.2e}")
open('r3_g13_out.txt', 'w').write('\n'.join(out) + '\n')
