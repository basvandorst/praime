# Independent constants: Artin S3 (d = -23, -31), Sym^2 E (11a, 37a, 15a), genus 2 (249, 277).
# Route 1 (direct): c = lim prod_{p<=P} F_p(1/p)/(1-1/p), F_p(X) = sum_k conj(lam(p^k)) mu(p^k) X^k from the local factor of L only;
#   log-Cesaro mean over [P/10, P], bias correction exp(nu/(sqrt p ln p)), nu = number of orthogonal irreducible pieces of Ad.
# Route 2 (exact): c = B(1)/L(1, Ad) with my own B_p = F_p(1/p) prod_{i,j} (1 - a_i conj(a_j)/p)^{-1} (all p, true local factors).
import os, sys, time, math, json, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6); pari.set_real_precision(30)
which = sys.argv[1:] or ['R23', 'R31', 'S11a', 'S37a', 'S15a', 'G249', 'G277']
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def local_F_and_B(poly):
    """poly: coefficients [1, b1, b2, ...] of the analytic local factor L_p(X)^{-1} = sum b_k X^k (complex ok).
    Returns F_p(X) coefficients (via lam(p^k) from 1/poly and mu(p^k) = b_k) and the Satake roots."""
    b = np.array(poly, dtype=complex); n = len(b) - 1
    lam = np.zeros(n + 1, complex); lam[0] = 1
    for k in range(1, n + 1):
        lam[k] = -sum(b[j] * lam[k - j] for j in range(1, k + 1))
    Fc = np.conj(lam) * b
    # roots alpha_i: poly(X) = prod (1 - alpha_i X)  -> alpha_i = 1/roots
    r = np.roots(b[::-1]) if n > 0 else np.array([])
    alpha = 1 / r
    return Fc, alpha

def Fval(Fc, X): return sum(Fc[k] * X ** k for k in range(len(Fc)))

def RSfactor(alpha, X):
    v = 1.0 + 0j
    for a in alpha:
        for b in alpha:
            v *= 1 / (1 - a * np.conj(b) * X)
    return v

def finish(key, pr, logterm, Bterm, nu, A1, extra):
    cum = np.cumsum(logterm)
    P = pr[-1]; m = pr >= P / 10; w = 1.0 / pr[m]
    d_end = math.exp(cum[-1]); d_sm = math.exp((cum[m] * w).sum() / w.sum())
    d_smb = math.exp(((cum[m] + nu / (np.sqrt(pr[m]) * np.log(pr[m]))) * w).sum() / w.sum())
    B = math.exp(np.sum(np.log(Bterm)))
    c_ex = B / A1 if A1 else None
    print('%-5s direct P=%d: end %.5f log-Cesaro %.5f +bias %.5f | B(1) = %.6f, L(1,Ad) = %s, exact c = %s %s'
          % (key, P, d_end, d_sm, d_smb, B, ('%.8f' % A1) if A1 else '-', ('%.6f' % c_ex) if c_ex else '-', extra))
    return dict(direct_end=d_end, direct_smooth=d_sm, direct_smooth_bias=d_smb, B=B, A1=A1, c_exact=c_ex)

res = {}
for key in which:
    t0 = time.time()
    if key in ('R23', 'R31'):
        pol = 'x^3-x-1' if key == 'R23' else 'x^3+x-1'; d = -23 if key == 'R23' else -31
        P = 10 ** 6; pr = primes(P)
        f = pari(pol)
        nroots = np.array([0 if p == -d else len(pari.polrootsmod(f, int(p))) for p in pr])
        lt = np.zeros(len(pr)); Bt = np.ones(len(pr))
        X = 1.0 / pr
        for i, p in enumerate(pr):
            x = X[i]
            if p == -d: lt[i] = 0.0; Bt[i] = 1 / (1 - x); continue
            nr = nroots[i]
            if nr == 3: lt[i] = math.log(1 - 3 * x); Bt[i] = (1 - 3 * x) / (1 - x) ** 3
            elif nr == 0: lt[i] = 0.0; Bt[i] = 1 / (1 - x ** 3)
            elif nr == 1: lt[i] = math.log(1 + x); Bt[i] = 1 / (1 - x * x)
            else: raise ValueError((p, nr))
        h = int(pari.qfbclassno(d)); bnf = pari.bnfinit(f, 1); hK = int(pari('bnfinit(%s,1).no' % pol)); RK = float(pari('bnfinit(%s,1).reg' % pol))
        A1 = (math.pi * h / math.sqrt(-d)) * (2 * math.pi * hK * RK / math.sqrt(-d))
        res[key] = finish(key, pr, lt, Bt, 2, A1, '(h(d)=%d, h_K=%d, R_K=%.6f; split/3-cycle/transp fractions %.3f %.3f %.3f)'
                          % (h, hK, RK, (nroots == 3).mean(), (nroots == 0).mean(), (nroots == 1).mean()))
    elif key.startswith('S'):
        ainv = {'S11a': [0, -1, 1, -10, -20], 'S37a': [0, 0, 1, -1, 0], 'S15a': [1, 1, 1, -10, -10]}[key]
        E = pari.ellinit(ainv); N = int(pari.ellglobalred(E)[0])
        P = 10 ** 6; pr = primes(P)
        an = pari.ellan(E, P)
        ap = np.array([int(an[int(p) - 1]) for p in pr], float)
        lt = np.zeros(len(pr)); Bt = np.ones(len(pr))
        for i, p in enumerate(pr):
            x = 1.0 / p
            if N % p == 0:
                # multiplicative: Sym^2 local factor (1 - X/p)^{-1}: F_p(X) = 1 - X/p^2 ; RS factor of St3 x St3 = zeta_p(s) zeta_p(s+1) zeta_p(s+2)
                Fv = 1 - x / p ** 2
                lt[i] = math.log(Fv / (1 - x)); Bt[i] = Fv / ((1 - x) * (1 - x / p) * (1 - x / p ** 2))
                continue
            lam = ap[i] / math.sqrt(p); s = lam * lam - 1
            Fc, alpha = local_F_and_B([1, -s, s, -1])
            Fv = Fval(Fc, x).real
            lt[i] = math.log(Fv / (1 - x)); Bt[i] = (Fv * RSfactor(alpha, x)).real
        # PARI L-values: analytic s = 1 <-> PARI s = 1 + m/2
        L2 = complex(pari.lfun(pari.lfunsympow(E, 2), 2)).real
        L4 = complex(pari.lfun(pari.lfunsympow(E, 4), 3)).real
        res[key] = finish(key, pr, lt, Bt, 2, L2 * L4, '(N=%d, L(1,Sym2)=%.6f, L(1,Sym4)=%.6f)' % (N, L2, L4))
    elif key.startswith('G'):
        PQ = {'G249': ('x^2+x', 'x^3+1'), 'G277': ('x^5-2*x^3+2*x^2-x', '1')}[key]
        Lg = pari('lfungenus2([%s,%s])' % PQ); N = int(pari('lfunparams(lfungenus2([%s,%s]))[1]' % PQ))
        P = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 30000
        pr = primes(P)
        lt = np.zeros(len(pr)); Bt = np.ones(len(pr))
        for i, p in enumerate(pr):
            p = int(p); x = 1.0 / p
            if N % p == 0:
                eu = pari.lfuneuler(Lg, p)          # 1/poly(T), arithmetic T = p^-s
                poly = pari.denominator(eu); num = pari.numerator(eu)
                cf = [float(pari.polcoef(poly, j)) / float(num) for j in range(int(pari.poldegree(poly)) + 1)]
                b = [cf[j] / p ** (j / 2) for j in range(len(cf))]
            else:
                cp = pari('hyperellcharpoly([Mod(1,%d)*(%s), Mod(1,%d)*(%s)])' % (p, PQ[0], p, PQ[1]))
                a1 = -int(pari.polcoef(cp, 3)); a2 = int(pari.polcoef(cp, 2))
                b = [1, -a1 / math.sqrt(p), a2 / p, -a1 / math.sqrt(p), 1]
            Fc, alpha = local_F_and_B(b)
            Fv = Fval(Fc, x).real
            lt[i] = math.log(Fv / (1 - x)); Bt[i] = (Fv * RSfactor(alpha, x)).real
            if N % p == 0: print('   bad p=%d analytic local poly %s, F_p(1/p)=%.6f' % (p, np.round(b, 5), Fv))
        res[key] = finish(key, pr, lt, Bt, 2, None, '(N=%d)' % N)
        res[key]['A1_implied_by_direct'] = res[key]['B'] / res[key]['direct_smooth_bias']
        print('      L(1,Ad) implied = B/c_direct = %.5f' % res[key]['A1_implied_by_direct'])
    print('      (%.0fs)' % (time.time() - t0)); sys.stdout.flush()
json.dump(res, open(OUT + 'r6b_%s.json' % '_'.join(k for k in which if not k.isdigit()), 'w'), indent=1)
