"""constants of the degree-4 follow-up (no zero computed). c by the direct product G = prod F_p(1/p)/(1 - 1/p): exact local
factors for p <= P0 (Sym^3, RS: 3000; genus 2: 1000), first order (1 - |lambda_p|^2/p)/(1 - 1/p) to P (1e5; genus 2 5e4),
log-Cesaro mean over P' in [P/10, P] with the prime-square bias nu/(sqrt P' ln P') (nu = 3 Sym^3, 2 genus 2, 3 RS). Sym^3 exact
route B/(L(1,Sym^2)L(1,Sym^4)L(1,Sym^6)) on the three smallest conductors as a check. RS: eps by the functional equation.
usage: c6_constants.py test|dry outdir"""
import sys, os, json, math, cmath, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np, cypari2
from bz4_lib import build, params
from ulib import to_c, local_poly, local_data, primes_upto, theta
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6); pari.default('realprecision', 38)
mode, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
def direct(lam, polyfn, P0, P, nu):
    pr = primes_upto(P); logs = []; Fp = {}
    for p in pr[pr <= P0]:
        p = int(p); Pp = polyfn(p); F, S, _ = local_data(Pp, p); Fp[p] = (F, S, Pp)
        logs.append(cmath.log(F) - math.log(1 - 1 / p))
    big = pr[pr > P0]
    lb = np.log(1 - np.abs(lam[big - 1]) ** 2 / big) - np.log(1 - 1 / big)
    cum = np.concatenate([np.cumsum(logs), sum(logs) + np.cumsum(lb)])
    sel = pr >= P // 10; w = 1.0 / pr[sel]; corr = nu / (np.sqrt(pr[sel]) * np.log(pr[sel]))
    c_s = complex(np.exp(np.sum((cum[sel] + corr) * w) / w.sum()))
    c_mid = complex(np.exp(cum[pr <= P // 3][-1]))
    m2 = float(np.mean(np.abs(lam[pr[pr > 100] - 1]) ** 2))
    return c_s.real, complex(np.exp(cum[-1])).real, c_mid.real, Fp, m2
def lines(eps, c, Fp):
    out = {}
    for p in (2, 3):
        F, S, Pp = Fp[p]
        for v in (1, -1): out['%s%d' % ('' if v > 0 else '1/', p)] = complex(-eps * c * S[v] / F)
    return out
def one(key, spec, T2, cls, exact_check=False):
    t0 = time.time()
    if spec[0] == 'rs':
        feqs = {}
        for e in (1, -1):
            feqs[e] = float(pari.lfuncheckfeq(build(pari, ['rs', spec[1], spec[2], e])))
        e_ok = min(feqs, key=lambda e: feqs[e]); spec = ['rs', spec[1], spec[2], e_ok]
    L = build(pari, spec); Q, k, mu = params(pari, L); eps = to_c(pari.lfunrootres(L)[2]); feq = float(pari.lfuncheckfeq(L))
    r = dict(key=key, cls=cls, spec=spec, Q=Q, k=k, mu=mu, degree=len(mu), eps=[eps.real, eps.imag], T2=T2, feq=feq,
             M=float(theta(T2, Q, mu)) / math.pi, M_band=float(theta(T2, Q, mu) - theta(T2 / 4, Q, mu)) / math.pi)
    if spec[0] == 'sym':
        an = pari.lfunan(L, 100000); lam = np.array([to_c(an[n - 1]) for n in range(1, 100001)]) / np.arange(1, 100001) ** ((k - 1) / 2)
        c_s, c_e, c_m, Fp, m2 = direct(lam, lambda p: local_poly(pari, L, p, k), 3000, 100000, 3)
        E = pari.ellinit(spec[1]); sgn = pari.lfunsympow(E, 3)
        Ls = pari('my(Lt = lfunsympow(ellinit(%s), 3)); Lt[6] = -Lt[6]; Lt' % spec[1])
        r['feq_flipped_eps'] = float(pari.lfuncheckfeq(Ls)); r['w_E'] = int(pari.ellrootno(E))
        if exact_check:
            objs = [pari.lfunsympow(E, 2 * j) for j in (1, 2, 3)]
            A1 = 1.0
            for j, O in enumerate(objs): A1 *= to_c(pari.lfun(O, 1 + (j + 1))).real
            lb = 0.0
            for p, (F, S, Pp) in Fp.items():
                a = np.array([1.0])
                for j, O in enumerate(objs): a = np.convolve(a, local_poly(pari, O, p, 2 * (j + 1) + 1))
                lb += math.log((F * (1 / sum(a[jj] * p ** (-jj) for jj in range(len(a))))).real / (1 - 1 / p))
            r['c_exact_check'] = math.exp(lb) / A1
    elif spec[0] == 'g2':
        an = pari.lfunan(L, 50000); lam = np.array([to_c(an[n - 1]) for n in range(1, 50001)]) / np.sqrt(np.arange(1, 50001))
        c_s, c_e, c_m, Fp, m2 = direct(lam, lambda p: local_poly(pari, L, p, k), 1000, 50000, 2)
    else:
        E1 = pari.ellinit(spec[1]); E2 = pari.ellinit(spec[2]); N1 = int(pari.ellglobalred(E1)[0]); N2 = int(pari.ellglobalred(E2)[0])
        a1 = pari.ellan(E1, 100000); a2 = pari.ellan(E2, 100000)
        ns = np.arange(1, 100001)
        lam = np.array([float(a1[n - 1]) * float(a2[n - 1]) for n in range(1, 100001)]) / ns     # correct at primes (only primes used)
        def poly(p):
            a = float(a1[p - 1]); b = float(a2[p - 1])
            if N1 % p == 0: P = [1, -a * b, a * a * p]
            elif N2 % p == 0: P = [1, -a * b, b * b * p]
            else: P = [1, -a * b, p * a * a + p * b * b - 2 * p * p, -a * b * p * p, p ** 4]
            return [c * p ** (-j) for j, c in enumerate(P)]
        c_s, c_e, c_m, Fp, m2 = direct(lam, poly, 3000, 100000, 3)
        r['feq_other_eps'] = feqs[-spec[3]]
    r.update(c=c_s, c_direct_end=c_e, c_direct_P3=c_m, mean_lam2=m2, c_method='direct (log-Cesaro + bias)',
             lines={kk: [v.real, v.imag] for kk, v in lines(eps, c_s, Fp).items()}, secs=round(time.time() - t0, 1))
    return r
items = []
if mode == 'test':
    TS = json.load(open('testset4.json'))
    sym = sorted(TS['sym3'], key=lambda v: v['N'])
    for i, v in enumerate(TS['sym3']):
        items.append(('S3_%s_N%d_%s' % (v['stratum'], v['N'], '_'.join(str(a).replace('-', 'm') for a in v['ainv'])), ['sym', v['ainv'], 3], v['T2'], 'sym3', v in sym[:3]))
    for i, g in enumerate(TS['g2']):
        items.append(('G2_N%d_%d' % (g[0], i), ['g2', g[1], g[2]], 150, 'g2', False))
    for i, (x, y) in enumerate(TS['rs']):
        items.append(('RS_%dx%d_%d' % (x[0], y[0], i), ['rs', x[1], y[1], 1], 60, 'rs', False))
else:
    R = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'blind', 'predictions_blind.json')))
    for key, e in R['L'].items():
        if e['cls'] == 'T4a': items.append((key, ['sym', e['spec'][1], 3], e['T2'], 'sym3', False))
        if e['cls'] == 'T4b': items.append((key, ['g2', e['spec'][1], e['spec'][2]], e['T2'], 'g2', False))
for key, spec, T2, cls, ex in items:
    fn = os.path.join(OUT, key + '.json')
    if os.path.exists(fn): continue
    r = one(key, spec, T2, cls, ex)
    json.dump(r, open(fn, 'w'), indent=1)
    print('%-34s c=%+.5f (end %+.5f, P/3 %+.5f%s) eps %+.0f feq %.0f%s M %.0f  %.0fs' % (key, r['c'], r['c_direct_end'], r['c_direct_P3'],
          (', exact %.5f' % r['c_exact_check']) if 'c_exact_check' in r else '', r['eps'][0], r['feq'],
          (' flipped %.0f w(E) %+d' % (r['feq_flipped_eps'], r['w_E'])) if 'feq_flipped_eps' in r else ((' other-eps %.0f' % r['feq_other_eps']) if 'feq_other_eps' in r else ''),
          r['M'], r['secs']), flush=True)
