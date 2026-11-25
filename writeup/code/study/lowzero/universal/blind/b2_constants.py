"""constants for the blind (or dry) set: c_L, eps, theta data, line factors, product-mechanism values. No zero is computed.
c: exact B(1)/L(1, Ad) when PARI gives L(1, Ad) with lfuncheckfeq <= -20 (mf: arithmetic GP Sym^2 with bad factors
1/(1 - a_p^2 X) at p | N; sym m: product of lfunsympow(E, 2j), j = 1..m), B over p <= 3000; always also the direct product
(exact local factors p <= P0, first order to P, log-Cesaro over [P/10, P], bias nu/(sqrt P ln P)); genus 2: direct only.
usage: b2_constants.py <testset.json> <outdir> [class ...]"""
import sys, os, json, math, cmath, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from bz_lib import build, params
from ulib import to_c, local_poly, local_data, series_inv, primes_upto, theta
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6); pari.default('realprecision', 38)
TS = json.load(open(sys.argv[1])); OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
CLS = sys.argv[3:]
NU = dict(mf=1, sym2=2, sym3=3, g2=2, artin=2)

def lam_vec(L, k, P):
    an = pari.lfunan(L, P)
    return np.array([to_c(an[n - 1]) for n in range(1, P + 1)]) / np.arange(1, P + 1) ** ((k - 1) / 2)

def direct(L, k, Q, polyfn, P0, P, nu, lam=None):
    lam = lam_vec(L, k, P) if lam is None else lam; pr = primes_upto(P)
    logs = []; Fp = {}
    for p in pr[pr <= P0]:
        p = int(p); Pp = polyfn(p); F, S, _ = local_data(Pp, p); Fp[p] = (F, S, Pp)
        logs.append(cmath.log(F) - math.log(1 - 1 / p))
    big = pr[pr > P0]
    lb = np.log(1 - np.abs(lam[big - 1]) ** 2 / big) - np.log(1 - 1 / big)
    cum = np.concatenate([np.cumsum(logs), sum(logs) + np.cumsum(lb)])
    sel = pr >= P // 10; w = 1.0 / pr[sel]
    corr = nu / (np.sqrt(pr[sel]) * np.log(pr[sel]))
    c_s = complex(np.exp(np.sum((cum[sel] + corr) * w) / w.sum()))
    m2 = float(np.mean(np.abs(lam[pr[pr > 100] - 1]) ** 2))
    return c_s, complex(np.exp(cum[-1])), Fp, m2

def exact(Fp, Ap_fn, A1, P0):
    lb = 0j
    for p, (F, S, Pp) in Fp.items():
        if p > P0: continue
        a = Ap_fn(p); Aval = 1 / sum(a[j] * p ** (-j) for j in range(len(a)))
        lb += cmath.log(F * Aval) - math.log(1 - 1 / p)
    return complex(cmath.exp(lb)) / A1

def lines(eps, c, Fp):
    out = {}
    for p in (2, 3):
        F, S, Pp = Fp[p]
        for v in (1, -1):
            out['%s%d' % ('' if v > 0 else '1/', p)] = complex(-eps * c * S[v] / F)
    return out

def one(key, spec, T2, cls):
    t0 = time.time()
    L = build(pari, spec); Q, k, mu = params(pari, L); eps = to_c(pari.lfunrootres(L)[2])
    r = dict(key=key, cls=cls, spec=spec, Q=Q, k=k, mu=mu, degree=len(mu), eps=[eps.real, eps.imag], T2=T2,
             M=float(theta(T2, Q, mu)) / math.pi)
    kind = spec[0]
    if kind == 'char':
        q = spec[1]; c = q / float(pari.eulerphi(q))
        Fp = {}
        for p in (2, 3):
            Pp = local_poly(pari, L, p, k); F, S, _ = local_data(Pp, p); Fp[p] = (F, S, Pp)
        r.update(c=c, c_method='q/phi(q) (theorem)')
    elif kind == 'mf':
        kk, N = spec[1], spec[2]
        pari('global_a = mfcoefs(mfeigenbasis(mfinit([%d,%d],0))[%d], 3000)' % (N, kk, spec[3] + 1))
        av = pari('global_a')
        lam = np.array([float(av[n]) / n ** ((kk - 1) / 2) for n in range(1, 3001)])
        def poly(p, lam=lam, N=N):
            return [1, -lam[p - 1]] if N % p == 0 else [1, -lam[p - 1], 1]
        c_s, c_e, Fp, m2 = direct(L, k, Q, poly, 3000, 3000, NU['mf'], lam=lam)
        pari('global_k = %d; global_N = %d' % (kk, N))
        try:
            S2 = pari('lfuncreate([n -> direuler(p=2, n, if(global_N%p, 1/((1-p^(global_k-1)*X)*(1-(global_a[p+1]^2-2*p^(global_k-1))*X+p^(2*global_k-2)*X^2)), 1/(1-global_a[p+1]^2*X))), 0, [2-global_k,0,1], 2*global_k-1, global_N^2, 1])')
            feq = float(pari.lfuncheckfeq(S2)); A1 = to_c(pari.lfun(S2, kk)).real
        except Exception as ex:
            feq = 99.0; A1 = float('nan')
        def Ap(p, lam=lam, N=N):
            if N % p == 0: return [1, -1.0 / p]
            l2 = abs(lam[p - 1]) ** 2; return list(np.convolve([1, -1], [1, -(l2 - 2), 1]))
        cx = exact(Fp, Ap, A1, 3000).real if feq <= -8 else float('nan')
        r.update(c_direct=c_s.real, c_direct_end=c_e.real, c_exact=cx, A_feq=feq, mean_lam2=m2)
        r.update(c=cx if feq <= -8 else c_s.real, c_method='exact (Sym2 GP, feq %.0f)' % feq if feq <= -8 else 'direct p<=3000 (Sym2 feq %.0f)' % feq)
    elif kind == 'sym':
        E = pari.ellinit(spec[1]); m = spec[2]
        objs = [pari.lfunsympow(E, 2 * j) for j in range(1, m + 1)]
        def poly(p, L=L, k=k): return local_poly(pari, L, p, k)
        c_s, c_e, Fp, m2 = direct(L, k, Q, poly, 3000, 100000, NU['sym%d' % m])
        def Ap(p):
            a = np.array([1.0])
            for j, O in enumerate(objs):
                a = np.convolve(a, local_poly(pari, O, p, 2 * (j + 1) + 1))
            return list(a)
        if m >= 3:
            r.update(c_direct=c_s.real, c_direct_end=c_e.real, mean_lam2=m2, c=c_s.real,
                     c_method='direct (exact p <= 3000, first order to 1e5, nu = 3; L(1, Sym^6 E) too costly)')
            r['lines'] = {kk2: [vv.real, vv.imag] for kk2, vv in lines(eps, r['c'], Fp).items()}
            r['secs'] = round(time.time() - t0, 1)
            return r
        A1 = 1.0; feqs = []
        for j, O in enumerate(objs):
            A1 *= to_c(pari.lfun(O, 1 + (j + 1))).real; feqs.append(float(pari.lfuncheckfeq(O)))
        cx = exact(Fp, Ap, A1, 3000).real
        r.update(c_direct=c_s.real, c_direct_end=c_e.real, c_exact=cx, A_feq=max(feqs), mean_lam2=m2)
        r.update(c=cx if max(feqs) <= -20 else c_s.real, c_method=('exact B/prod L(1,Sym^2j E), feq %.0f' % max(feqs)) if max(feqs) <= -20 else 'direct')
    elif kind == 'g2':
        def poly(p, L=L, k=k): return local_poly(pari, L, p, k)
        c_s, c_e, Fp, m2 = direct(L, k, Q, poly, 1000, 50000, NU['g2'])
        r.update(c_direct=c_s.real, c_direct_end=c_e.real, mean_lam2=m2, c=c_s.real,
                 c_method='direct (exact p <= 1000, first order to 5e4, nu = 2)')
    elif kind == 'nf':
        r.update(c=0.0, c_method='0 (product zeta * L(rho))')
        Lr = build(pari, ['artin', spec[1]]); Qr, kr, mur = params(pari, Lr)
        d = int(pari.nfdisc(pari(spec[1]))); Lc = pari.lfuncreate(d)
        def poly(p, Lr=Lr): return local_poly(pari, Lr, p, 1)
        c_s, c_e, Fp, m2 = direct(Lr, 1.0, Qr, poly, 3000, 100000, NU['artin'])
        def Ap(p): return list(np.convolve(local_poly(pari, Lc, p, 1), local_poly(pari, Lr, p, 1)))
        A1 = to_c(pari.lfun(Lc, 1)).real * to_c(pari.lfun(Lr, 1)).real
        cr = exact(Fp, Ap, A1, 3000).real
        epr = to_c(pari.lfunrootres(Lr)[2])
        K = int(T2 / (2 * math.pi))
        # mechanism: rho(1/k) multiplicative, psi_rho
        rho = {}; psi = 0.0; psiH = 0.0; lnk = {}
        for p in primes_upto(K):
            p = int(p); Pp = local_poly(pari, Lr, p, 1); lam = series_inv(Pp, 80)
            F = sum(np.conj(lam[j]) * Pp[j] * p ** (-j) for j in range(len(Pp)))
            al = [1 / x for x in np.roots(Pp[::-1])] if len(Pp) > 1 else []
            v = 1
            while p ** v <= K:
                h = sum(np.conj(lam[hh - v]) * Pp[hh] * p ** (-(2 * hh - v) / 2) for hh in range(v, len(Pp)))
                rho[(p, v)] = p ** (v / 2) * h / F
                lv = sum(a ** v for a in al).real * math.log(p)
                psi += lv; psiH += lv * math.sin(math.pi * 2 * math.pi * p ** v / T2) ** 2
                v += 1
        Msum = 0; MsumH = 0
        for kk in range(1, K + 1):
            f = pari.factor(kk); val = 1
            for i in range(len(f[0])): val *= rho[(int(f[0][i]), int(f[1][i]))]
            Msum += val; MsumH += val * math.sin(math.pi * 2 * math.pi * kk / T2) ** 2
        r.update(d=d, rho=dict(Q=Qr, mu=mur, eps=[epr.real, epr.imag], c=cr, c_direct=c_s.real, mean_lam2=m2,
                               M=float(theta(T2, Qr, mur)) / math.pi, lines={kk2: [vv.real, vv.imag] for kk2, vv in lines(epr, cr, Fp).items()}),
                 mech=dict(K=K, sum_rho_1_over_k=complex(Msum).real, sum_rho_1_over_k_hann=complex(MsumH).real, psi_rho=psi, psi_rho_hann=psiH,
                           M1=-cr * complex(Msum).real, M1_hann=-cr * complex(MsumH).real, M2=-psi, M2_hann=-psiH, M3=0.0))
        r['secs'] = round(time.time() - t0, 1)
        return r
    if kind != 'nf':
        r['lines'] = {kk2: [vv.real, vv.imag] for kk2, vv in lines(eps, r['c'], Fp).items()}
    r['secs'] = round(time.time() - t0, 1)
    return r

def items(cls):
    out = []
    for x in TS[cls]['test']:
        if cls == 'T1':
            out += [('T1_%d_%d' % (x[0], x[1]), ['char', x[0], x[1]], 1000), ('T1_%d_%d' % (x[0], x[2]), ['char', x[0], x[2]], 1000)]
        elif cls == 'T2': out.append(('T2_k%d_N%d_%d' % (x[0], x[1], x[2]), ['mf', x[0], x[1], x[2]], TS.get('height', {}).get('T2', 500)))
        elif cls == 'T3': out.append(('T3_N%d_%s' % (x[0], '_'.join(str(a).replace('-', 'm') for a in x[1])), ['sym', x[1], 2], TS.get('height', {}).get('T3', 200)))
        elif cls == 'T4a': out.append(('T4a_N%d_%s' % (x[0], '_'.join(str(a).replace('-', 'm') for a in x[1])), ['sym', x[1], 3], TS.get('height', {}).get('T4a', 100)))
        elif cls == 'T4b': out.append(('T4b_N%d_%d' % (x[0], TS['T4b']['test'].index(x)), ['g2', x[1], x[2]], TS.get('height', {}).get('T4b', 150)))
        elif cls == 'T5': out.append(('T5_d%d' % (-x[0]), ['nf', x[1]], 300))
    return out
for cls in CLS:
    for key, spec, T2 in items(cls):
        fn = os.path.join(OUT, key + '.json')
        if os.path.exists(fn): continue
        try:
            r = one(key, spec, T2, cls)
            json.dump(r, open(fn, 'w'), indent=1)
            print('%-28s c=%+.5f [%s] eps=%+.3f%+.3fi M=%.0f %s %.0fs' % (key, r['c'] if 'c' in r else 0, r.get('c_method', '')[:34], r['eps'][0], r['eps'][1], r['M'],
                  ('c_rho=%.4f M1=%.1f M2=%.1f' % (r['rho']['c'], r['mech']['M1'], r['mech']['M2'])) if 'rho' in r else
                  ('direct %.5f' % r.get('c_direct', float('nan'))), r['secs']), flush=True)
        except Exception as e:
            import traceback; traceback.print_exc(); print(key, 'FAILED', str(e)[:200], flush=True)
