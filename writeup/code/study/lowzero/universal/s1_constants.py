"""constants of the universal one-point law, computed from Euler products and L-values only (no zero used).
c_L = lim_{x->0+} F(1+x)/x, F(s) = sum conj(lambda(n)) mu_L(n) n^-s.
 (a) direct: prod_p F_p(1/p)/(1-1/p) (conditionally convergent), exact local factors p <= P0, first order |lambda_p|^2 to P,
     log-Cesaro mean over P' in [P/10, P] with the prime-square bias correction nu/(sqrt(P') ln P').
 (b) exact: prod_{p<=P0} F_p(1/p) A_p(1/p)/(1-1/p) / A(1), A = L(s, Ad) (absolutely convergent product; A(1) from PARI).
Lines: amplitude per unit (T2-T1)/2pi of sum e^{2i theta(gamma)} r^{-i gamma}, r = p^v:  -eps c S_p(v)/F_p(1/p).
usage: s1_constants.py [keys...]  -> constants_<key>.json"""
import sys, os, time, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from ulib import *
pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6); pari.default('realprecision', 38)
P0, P = 3000, 10 ** 5
NU = {'GL2 level 1': 1, 'GL2 weight 2 irrational': 1, 'GL2 nebentypus (non-self-dual)': 1, 'GL3 Sym2 E': 2,
      'GL2 Artin S3': 2, 'GL4 genus 2': 2}

def mult_coeffs(M, locfn):
    """coefficients of prod_p 1/locfn(p)(p^-s) up to M (analytic); locfn(p) -> polynomial coefficient list."""
    c = np.zeros(M + 1, complex); c[1] = 1
    spf = np.zeros(M + 1, np.int64)
    for p in primes_upto(M):
        m = spf[p::p]; m[m == 0] = p; spf[p::p] = m
    loc = {}
    for n in range(2, M + 1):
        p = int(spf[n]); m = n; j = 0
        while m % p == 0: m //= p; j += 1
        if p not in loc:
            J = int(math.log(M) / math.log(p)) + 1
            loc[p] = series_inv(locfn(p), J)
        c[n] = loc[p][j] * c[m]
    return c[1:]

def make_lfun(coeffs, mu, Q, real=True):
    v = pari([float(x.real) for x in coeffs]) if real else pari([complex(x) for x in coeffs])
    return pari.lfuncreate([v, 0 if real else 1, [float(m) for m in mu], 1, Q, 1])

def adjoint_object(key, spec, Q, k, loc):
    """A = L(s, Ad pi) with local polynomial function (analytic) and value A(1); None if not available."""
    fam = spec[0]
    if fam == 'mf':
        N = spec[1]
        def Ap(p):
            if N % p == 0:
                return [1, -1.0 / p] if spec[3] == 1 else [1, -1.0]          # Steinberg: 1/(1-X/p); ramified principal series: 1/(1-X)
            l2 = abs(loc.lam[p - 1]) ** 2
            return list(np.convolve([1, -1], [1, -(l2 - 2), 1]))
        M = 40000
        co = mult_coeffs(M, Ap)
        A = make_lfun(co, [1, k - 1, k], N * N)
        feq = float(pari.lfuncheckfeq(A))
        val = to_c(pari.lfun(A, 1))
        return Ap, val.real, feq, 'constructed Sym2/Ad (analytic, conductor %d, mu = [1,%g,%g], lfuncheckfeq %.0f)' % (N * N, k - 1, k, feq)
    if fam == 'sym':
        E = pari.ellinit(spec[1]); S2 = pari.lfunsympow(E, 2); S4 = pari.lfunsympow(E, 4)
        def Ap(p):
            return list(np.convolve(local_poly(pari, S2, p, 3), local_poly(pari, S4, p, 5)))
        val = to_c(pari.lfun(S2, 2)) * to_c(pari.lfun(S4, 3))
        return Ap, val.real, float(pari.lfuncheckfeq(S4)), 'L(1,Sym2 E) L(1,Sym4 E) (PARI lfunsympow)'
    if fam == 'artin':
        d = int(pari.nfdisc(pari(spec[1])))
        Lc = pari.lfuncreate(d); Lr = build(pari, spec)
        def Ap(p):
            return list(np.convolve(local_poly(pari, Lc, p, 1), local_poly(pari, Lr, p, 1)))
        val = to_c(pari.lfun(Lc, 1)) * to_c(pari.lfun(Lr, 1))
        return Ap, val.real, 0.0, 'L(1,chi_d) L(1,rho), d = %d' % d
    return None

def run(key, spec, T2, fam):
    t0 = time.time()
    L = build(pari, spec)
    Q, k, mu = params(pari, L)
    eps = to_c(pari.lfunrootres(L)[2])
    res = dict(key=key, family=fam, spec=[str(x) for x in spec], Q=Q, k=k, mu=mu, degree=len(mu), eps=[eps.real, eps.imag], T2=T2)
    if spec[0] == 'dir':
        q = spec[1]; c = q / float(pari.eulerphi(q))
        res.update(c=c, c_method='exact q/phi(q) (degree-1 theorem; B = prod_{p|q}(1-1/p)^-1, Ad trivial)')
    elif fam.startswith('product'):
        res.update(c=0.0, c_method='0: L is a product, F has a zero of order >= 2 at s = 1')
    loc = Local(pari, L, spec, Q, k, P)
    pr = primes_upto(P)
    lam_p = loc.lam[pr - 1]
    res['mean_abs_lambda_p_sq'] = float(np.mean(np.abs(lam_p) ** 2))
    # exact local factors p <= P0
    logs = []; Fp = {}; Sp = {}; polys = {}
    for p in pr[pr <= P0]:
        p = int(p); Pp = loc.poly(p); polys[p] = Pp
        F, S, lam = local_data(Pp, p)
        Fp[p] = F; Sp[p] = S
        logs.append(cmath.log(F) - math.log(1 - 1 / p))
    logs = np.array(logs)
    # first order beyond P0
    big = pr[pr > P0]
    lb = np.log(1 - np.abs(loc.lam[big - 1]) ** 2 / big) - np.log(1 - 1 / big)
    cum = np.concatenate([np.cumsum(logs), logs.sum() + np.cumsum(lb)])
    nu = NU.get(fam, 0)
    sel = pr >= P // 10
    wts = 1.0 / pr[sel]
    corr = nu / (np.sqrt(pr[sel]) * np.log(pr[sel]))
    c_dir_end = complex(np.exp(cum[-1]))
    c_dir_smooth = complex(np.exp(np.sum((cum[sel] + corr) * wts) / wts.sum()))
    c_dir_1e4 = complex(np.exp(cum[pr <= 10 ** 4][-1]))
    res.update(c_direct_end=[c_dir_end.real, c_dir_end.imag], c_direct_smooth=[c_dir_smooth.real, c_dir_smooth.imag],
               c_direct_1e4=[c_dir_1e4.real, c_dir_1e4.imag], nu=nu)
    adj = None
    if 'c' not in res:
        adj = adjoint_object(key, spec, Q, k, loc)
    if adj is not None:
        Ap, A1, feqA, how = adj
        lb2 = []
        for p in pr[pr <= P0]:
            p = int(p); a = Ap(p)
            Aval = 1 / sum(a[j] * p ** (-j) for j in range(len(a)))
            lb2.append(cmath.log(Fp[p] * Aval) - math.log(1 - 1 / p))
        lb2 = np.array(lb2)
        Bc = complex(np.exp(lb2.sum()))
        c_ex = Bc / A1
        c_ex_1000 = complex(np.exp(lb2[pr[pr <= P0] <= 1000].sum())) / A1
        res.update(c_exact=[c_ex.real, c_ex.imag], c_exact_P0_1000=[c_ex_1000.real, c_ex_1000.imag], A1=A1, B=[Bc.real, Bc.imag],
                   A_how=how, A_feq=feqA)
        res.update(c=c_ex.real, c_method='exact: B(1)/L(1,Ad), B over p <= %d' % P0)
    elif 'c' not in res:
        res.update(c=c_dir_smooth.real, c_method='direct product, exact p <= %d, first order to %d, log-Cesaro over [P/10,P], bias nu = %d' % (P0, P, nu))
    c = res['c']
    lines = {}
    for p in (2, 3):
        for v in (1, -1):
            amp = -eps * c * Sp[p][v] / Fp[p]
            lines['%s%d' % ('' if v > 0 else '1/', p)] = [amp.real, amp.imag]
    res['lines_per_unit'] = lines
    res['F_p'] = {str(p): [Fp[p].real, Fp[p].imag] for p in (2, 3, 5, 7)}
    res['poly_2_3'] = {str(p): [[x.real, x.imag] for x in polys[p]] for p in (2, 3)}
    res['secs'] = round(time.time() - t0, 1)
    json.dump(res, open('constants_%s.json' % key, 'w'), indent=1)
    print('%-5s c=%.5f  [%s]  direct: end %.5f smooth %.5f 1e4 %.5f  exact %s  mean|lam_p|^2 %.4f  eps %.4f%+.4fi  (%.0fs)' % (
        key, c, res['c_method'][:20], c_dir_end.real, c_dir_smooth.real, c_dir_1e4.real,
        ('%.5f (P0=1000: %.5f, A-feq %.0f)' % (res['c_exact'][0], res['c_exact_P0_1000'][0], res['A_feq'])) if 'c_exact' in res else '-',
        res['mean_abs_lambda_p_sq'], eps.real, eps.imag, res['secs']), flush=True)

keys = sys.argv[1:]
for key, spec, T2, fam in CATALOGUE:
    if keys and key not in keys: continue
    try:
        run(key, spec, T2, fam)
    except Exception as e:
        import traceback; traceback.print_exc(); print(key, 'FAILED', str(e)[:300], flush=True)
