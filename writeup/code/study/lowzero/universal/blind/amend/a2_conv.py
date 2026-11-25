# Amended product mechanism for zeta_K = zeta * L(rho) (complex cubic fields), T2 = 300.
# rho-zeros weighted by e^{2i theta_K} = e^{2i theta_rho} e^{2i theta_zeta}: every line r = a/b of the degree-2 law of L(rho),
#   amplitude A_r = -eps c rho(r)/sqrt(ab) per dt/2pi times r^{it}, is swept by the zeta chirp e^{2i theta_zeta(t)};
#   contribution = A_r * (1/2pi) int_0^T e^{2i theta_zeta(t)} r^{it} dt (exact integral; stationary point t = 2 pi b/a, phase e(-b/a)).
# The registered (round-15x) form keeps only a = 1 in the stationary-phase approximation: -c sum_{k <= T/2pi} rho(1/k).
# Reads only predictions_blind.json (registered constants/predictions); never the zeros folder.
import os, sys, math, json, hashlib, time, numpy as np, cypari2
from scipy.special import loggamma
from math import gcd
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
BL = os.path.dirname(os.path.abspath(__file__))+'/../../../../../../../writeup/code/study/lowzero/universal/blind/'
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
PB = json.load(open(BL + 'predictions_blind.json'))
print('predictions_blind.json sha256', hashlib.sha256(open(BL + 'predictions_blind.json', 'rb').read()).hexdigest())
T = 300.0; x = T / (2 * math.pi); AMAX = 12; ATAIL = 30; PMAX = 50000
t0 = time.time()
# ---------------------------------------------------------------- field-independent integrals
tt = np.linspace(0.0, T, 300001); tt[0] = 1e-9; dt = tt[1] - tt[0]
thz = 0.5 * tt * math.log(1 / math.pi) + np.imag(loggamma((0.5 + 1j * tt) / 2))
E = np.exp(2j * thz); Wh = np.sin(math.pi * tt / T) ** 2
rats = []
for a in range(1, ATAIL + 1):
    for b in range(1, int(1.25 * a * x) + 4):
        if gcd(a, b) == 1: rats.append((a, b))
I = {}; Ih = {}
for (a, b) in rats:
    f = E * np.exp(1j * tt * math.log(a / b))
    I[(a, b)] = ((f.sum() - 0.5 * (f[0] + f[-1])) * dt) / (2 * math.pi)
    fh = f * Wh
    Ih[(a, b)] = ((fh.sum() - 0.5 * (fh[0] + fh[-1])) * dt) / (2 * math.pi)
print('%d integrals (a <= %d) in %.0fs' % (len(rats), ATAIL, time.time() - t0)); sys.stdout.flush()
# ---------------------------------------------------------------- per field local data
s = np.ones(PMAX + 1, bool); s[:2] = False
for i in range(2, int(PMAX ** .5) + 1):
    if s[i]: s[i * i::i] = False
primes = np.nonzero(s)[0]
def fac(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def series_inv(P, n):
    lam = [1.0] + [0.0] * n
    for k in range(1, n + 1):
        lam[k] = -sum(P[j] * lam[k - j] for j in range(1, min(k, len(P) - 1) + 1))
    return lam
def field(pol, d):
    nf = pari.nfinit(pari(pol)); D = int(pari.nfdisc(pari(pol))); assert D == d, (D, d)
    pd = abs(int(pari.poldisc(pari(pol))))
    loc = {}
    def local(p):
        if p in loc: return loc[p]
        if pd % p == 0:
            dec = pari.idealprimedec(nf, p)
            Pz = np.poly1d([1.0])
            for pr in dec:
                f_ = int(pr[3]); Pz = Pz * np.poly1d([-1.0] + [0.0] * (f_ - 1) + [1.0])   # (1 - X^f) as poly in X (highest first)
            # divide by (1 - X)
            q, r = np.polydiv(Pz.coeffs, np.array([-1.0, 1.0]))
            assert np.allclose(r, 0)
            P = list(q[::-1])
        else:
            nr = len(pari.polrootsmod(pari(pol), p))
            P = {3: [1.0, -2.0, 1.0], 0: [1.0, 1.0, 1.0], 1: [1.0, 0.0, -1.0]}[nr]
        while len(P) > 1 and abs(P[-1]) < 1e-15: P = P[:-1]
        lam = series_inv(P, 80)
        mu = lambda h: P[h] if h < len(P) else 0.0
        X = 1.0 / p
        Fp = sum(lam[k] * mu(k) * X ** k for k in range(len(P)))
        chi = int(pari.kronecker(d, p))
        Lrr = (1 / (1 - X)) * (1 / (1 - chi * X)) * series_eval(lam, X)
        Bp = Fp * Lrr
        def S(v):
            if v >= 0: return sum(lam[v + h] * mu(h) * X ** h for h in range(len(P)))
            return sum(lam[m] * mu(-v + m) * X ** m for m in range(0, 70) if -v + m < len(P))
        loc[p] = dict(P=P, Fp=Fp, Lrr=Lrr, Bp=Bp, S=S)
        return loc[p]
    def series_eval(lam, X): return sum(lam[k] * X ** k for k in range(len(lam)))
    globals()['series_eval'] = series_eval
    logB = 0.0; zeros = []
    for p in primes:
        L_ = local(int(p))
        if abs(L_['Fp']) < 1e-12: zeros.append(int(p)); continue
        logB += math.log(L_['Bp']) if L_['Bp'] > 0 else float('nan')
    h = int(pari.qfbclassno(d)); bnf = pari.bnfinit(pari(pol), 1)
    hK = int(pari('bnfinit(%s,1).no' % pol)); RK = float(pari('bnfinit(%s,1).reg' % pol))
    A1 = (math.pi * h / math.sqrt(-d)) * (2 * math.pi * hK * RK / math.sqrt(-d))
    Breg = math.exp(logB)
    def crho(a, b):
        fa = fac(a); fb = fac(b); ps = set(fa) | set(fb)
        if any(z not in ps for z in zeros): return 0.0
        v = Breg / A1
        for p in ps:
            L_ = local(p); e_ = fa.get(p, 0) - fb.get(p, 0)
            v *= L_['S'](e_) * L_['Lrr'] / (L_['Bp'] if p not in zeros else 1.0)
        return v
    c = 0.0 if zeros else Breg / A1
    return dict(c=c, crho=crho, zeros=zeros, h=h, hK=hK, RK=RK, A1=A1)
globals()['series_eval'] = lambda lam, X: sum(lam[k] * X ** k for k in range(len(lam)))

def run(key, pol, d, eps, reg):
    F = field(pol, d)
    lines = {r: -eps * F['crho'](*r) / math.sqrt(r[0] * r[1]) for r in rats}
    amp = {}
    tot = 0j; toth = 0j; bya = {}; tail = 0j
    for r in rats:
        v = lines[r] * I[r]; vh = lines[r] * Ih[r]
        if r[0] <= AMAX: tot += v; toth += vh; bya[r[0]] = bya.get(r[0], 0j) + v
        else: tail += v
    simple = -eps * sum(F['crho'](1, k) for k in range(1, int(x) + 1))     # registered a = 1 stationary form
    chk = {kk: -eps * F['crho'](*rr) / math.sqrt(rr[0] * rr[1]) for kk, rr in (('2', (2, 1)), ('1/2', (1, 2)), ('3', (3, 1)), ('1/3', (1, 3)))}
    pa = {}
    for rr in rats: pa[rr[0]] = pa.get(rr[0], 0j) + lines[rr] * I[rr]
    partial = []; acc = 0j
    for A in range(1, ATAIL + 1):
        acc += pa.get(A, 0j)
        if A in (1, 2, 3, 6, 12, 18, 24, 30): partial.append((A, acc))
    return dict(partial=[(A, v) for A, v in partial], c=F['c'], zeros=F['zeros'], M1rho=tot, M1rho_hann=toth, bya=bya, tail=tail, simple=simple, linechk=chk, h=F['h'], hK=F['hK'], RK=F['RK'])

res = {}
# controls from round 1 (not in the blind set): must reproduce -44.58-8.33i and -22.26-0.93i
for key, pol, d in (('R23', 'x^3-x-1', -23), ('R31', 'x^3+x-1', -31)):
    r = run(key, pol, d, 1.0, None)
    print('control %s: c = %.6f, rho-zeros/theta_K (a<=12) = %.2f%+.2fi, a=1 simple %.2f' % (key, r['c'], r['M1rho'].real, r['M1rho'].imag, r['simple']))
for key, e in PB['products'].items():
    if e['cls'] != 'T5': continue
    pol = e['spec'][1]; d = int(e['d']); eps = complex(*e['rho']['eps']).real
    r = run(key, pol, d, eps, e)
    rd = e['readings']
    reg = rd['M1']['pred']
    lc = e['rho']['lines']
    dl = max(abs(r['linechk'][k] - complex(*lc[k])) for k in lc)
    res[key] = dict(d=d, pol=pol, c_rho_registered=e['rho']['c'], c_rho_own=r['c'], zeros_of_F=r['zeros'], line_check_maxdiff=dl,
                    M1_registered=reg, M1_registered_own_recompute=r['simple'],
                    M1_rho_part_amended=[r['M1rho'].real, r['M1rho'].imag], M3_amended=0.0,
                    M1_amended=[r['M1rho'].real, r['M1rho'].imag], M1_amended_hann=[r['M1rho_hann'].real, r['M1rho_hann'].imag],
                    by_a={str(a): [v.real, v.imag] for a, v in sorted(r['bya'].items())}, tail_a13_18=[r['tail'].real, r['tail'].imag],
                    sd_M1=rd['M1']['sd'], sd_M3=rd['M3']['sd'], M_total=e['M'], h_d=r['h'], h_K=r['hK'], R_K=r['RK'])
    ps = []; acc = 0j
    allb = {}
    for rr in rats:
        allb[rr[0]] = allb.get(rr[0], 0j) + lines[rr] * I[rr] if False else allb.get(rr[0], 0j)
    print('%-8s' % key, end=' ')
    res[key]['partial'] = [[A, v.real, v.imag] for A, v in r['partial']]
    print(' '.join('A=%d:%+.1f' % (A, v.real) for A, v in r['partial']))
    sys.stdout.flush()
json.dump(res, open(OUT + 'a2_conv.json', 'w'), indent=1)
print('done (%.0fs)' % (time.time() - t0))
