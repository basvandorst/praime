"""independent recomputation (no import of ulib or s4_sums):
(1) predictions.json hash; every zero file written after the freeze time;
(2) theta by mpmath.loggamma, the sharp sums and Hann sums by plain Python, against results.json;
(3) zero counts against theta(T)/pi + offset (the offset 1/2 for eps = +1, 0 for eps = -1, 1 for zeta-type poles);
(4) constants by other routes: Delta and the level-1 forms via the Petersson norm, L(1,Sym^2 f) = (pi/2)(4pi)^k <f,f>/Gamma(k), and
    B by an independent loop with PARI's ramanujantau/mfcoefs; R23/R31 via the class-number closed form; S11a by the direct product
    with ellap to 3*10^5; Dirichlet q/phi(q).
usage: validate.py <sha256>"""
import sys, os, json, math, cmath, hashlib, glob
import numpy as np, mpmath, cypari2
mpmath.mp.dps = 30
raw = open('predictions.json').read(); h = hashlib.sha256(raw.encode()).hexdigest()
print('[1] sha256 matches:', h == sys.argv[1], h)
P = json.loads(raw)
tfreeze = os.path.getmtime('predictions.json')
late = [f for f in glob.glob('zeros/zeros_*.npz') if os.path.getmtime(f) < tfreeze]
print('[1] zero files older than predictions.json:', late if late else 'none', '; freeze file mtime', tfreeze)
R = json.load(open('results.json'))
def th(t, Q, mu):
    s = mpmath.mpf(t) / 2 * mpmath.log(mpmath.mpf(Q) / mpmath.pi ** len(mu))
    for m in mu:
        s += mpmath.im(mpmath.loggamma((mpmath.mpf('0.5') + m + 1j * mpmath.mpf(t)) / 2))
    return s
worst = 0; worsth = 0
SKIP2 = len(sys.argv) > 2 and sys.argv[2] == 'skip2'
print('[2] sums recomputed with mpmath theta (|diff| sharp, Hann):')
for key, r in ([] if SKIP2 else R['per_L'].items()):
    e = P['Lfunctions'][key]; g = np.load('zeros/zeros_%s.npz' % key)['g']; T2 = e['T2']
    g = [float(x) for x in g if 0 < x <= T2]
    S = 0j; Sh = 0j
    for x in g:
        z = complex(mpmath.exp(2j * th(x, e['Q'], e['mu'])))
        S += z; Sh += math.sin(math.pi * x / T2) ** 2 * z
    d1 = abs(S - complex(*r['G'])); d2 = abs(Sh - complex(*r['G_hann'])); worst = max(worst, d1); worsth = max(worsth, d2)
    M = float(th(T2, e['Q'], e['mu'])) / math.pi
    print('   %-5s n=%d  M(T2)=theta/pi=%.2f  n - M = %+.2f   diff %.1e %.1e' % (key, len(g), M, len(g) - M, d1, d2))
print('   worst', worst, worsth)
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6); pari.default('realprecision', 38)
def primes(n): return [int(p) for p in pari.primes([2, n])]
print('[4] constants by other routes:')
for key, k in ([] if SKIP2 else (('D12', 12), ('F16', 16), ('F18', 18), ('F20', 20), ('F22', 22), ('F26', 26))):
    mf = pari('mfinit([1,%d],0)' % k); F = pari.mfeigenbasis(mf)[0]
    pet = float(pari.mfpetersson(pari.mfsymbol(mf, F)))
    LS = math.pi / 2 * (4 * math.pi) ** k * pet / math.gamma(k)
    co = pari.mfcoefs(F, 3000)
    lnB = 0.0
    for p in primes(3000):
        l = float(co[p]) / p ** ((k - 1) / 2); e1 = l * l - 1
        lnB += math.log((1 - e1 / p) / (1 - e1 / p + e1 / p ** 2 - 1 / p ** 3))
    c = math.exp(lnB) / LS
    print('   %-5s c (Petersson route, B over p<=3000) = %.6f   frozen %.6f' % (key, c, P['Lfunctions'][key]['c']))
for key, pol in (('R23', 'x^3-x-1'), ('R31', 'x^3+x-1')):
    hK = int(pari('bnfinit(%s,1).no' % pol)); RK = float(pari('bnfinit(%s,1).reg' % pol)); d = int(pari.nfdisc(pari(pol)))
    hd = int(pari.qfbclassno(d))
    LAd = (math.pi * hd / math.sqrt(-d)) * (2 * math.pi * hK * RK / math.sqrt(-d))
    # B: local factors from splitting of p in K (independent of the lfun objects)
    lnB = 0.0
    for p in primes(3000):
        dec = pari.idealprimedec(pari.nfinit(pari(pol)), p) if p < 50 or d % p == 0 else None
        # ζ_K Euler factor -> rho factor = zeta_K factor / zeta factor
        fs = [int(pr[3]) for pr in pari.idealprimedec(pari.nfinit(pari(pol)), p)] if (d % p == 0 or p < 0) else None
        if d % p == 0:
            Pz = np.array([1.0])
            for f in fs: Pz = np.convolve(Pz, np.r_[1.0, np.zeros(f - 1), -1.0])
            Pr = np.polydiv(Pz[::-1], np.array([1.0, -1.0])[::-1])[0][::-1]     # divide by (1 - X)
        else:
            nroots = len(pari.polrootsmod(pari(pol), p))
            Pr = {3: np.array([1.0, -2.0, 1.0]), 1: np.array([1.0, 0.0, -1.0]), 0: np.array([1.0, 1.0, 1.0])}[nroots]
        lam = [1.0]
        for j in range(1, 40):
            lam.append(-sum(Pr[i] * lam[j - i] for i in range(1, min(j, len(Pr) - 1) + 1)))
        X = 1.0 / p
        F = sum(lam[j] * Pr[j] * X ** j for j in range(len(Pr)))
        chi = int(pari.kronecker(d, p))
        A = 1 / ((1 - chi * X) * sum(Pr[j] * X ** j for j in range(len(Pr))))     # L_p(chi_d) L_p(rho)
        lnB += math.log(F * A / (1 - X))
    c = math.exp(lnB) / LAd
    print('   %-5s h(d)=%d hK=%d RK=%.6f  L(1,Ad)=%.6f  c (closed form) = %.6f   frozen %.6f' % (key, hd, hK, RK, LAd, c, P['Lfunctions'][key]['c']))
E = pari.ellinit([0, -1, 1, -10, -20]); lnF = 0.0; part = []
ps = primes(300000)
for p in ps:
    if p == 11:
        lnF += math.log((1 - 1 / p ** 3) / (1 - 1 / p)); continue      # Sym2 at multiplicative p: L_p = 1/(1 - X/p): lambda(p) = 1/p, mu(p) = -1/p, F_p(X) = 1 - X/p^2, at X = 1/p: 1 - 1/p^3
    ap = int(pari.ellap(E, p)); l2 = ap * ap / p
    a = l2 - 1                                                          # Satake of Sym2: {alpha^2, 1, alpha^-2}, lambda(p) = a
    P3 = [1.0, -a, a, -1.0]                                             # (1 - alpha^2 X)(1 - X)(1 - alpha^-2 X)
    lam = [1.0]
    for j in range(1, 30):
        lam.append(-sum(P3[i] * lam[j - i] for i in range(1, min(j, 3) + 1)))
    X = 1.0 / p
    F = sum(lam[j] * P3[j] * X ** j for j in range(4))
    lnF += math.log(F / (1 - X))
    if p > 30000: part.append(lnF)
print('   S11a  direct product to 3e5: end %.5f, log-mean over [3e4,3e5] %.5f   frozen (exact B/L(1,Ad)) %.5f' % (math.exp(lnF), math.exp(np.mean(part)), P['Lfunctions']['S11a']['c']))
for key, q in (('X7', 7), ('X5', 5)):
    print('   %-5s q/phi(q) = %.6f  frozen %.6f' % (key, q / float(pari.eulerphi(q)), P['Lfunctions'][key]['c']))
