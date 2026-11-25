"""Frozen predictions for the swap2 re-test. Reads ONLY the keys of the reserved zero files (no zero value is used).
Writes pred_OO.npz (546 real L-functions, all 148,785 pairs; HI = r134_C uses the R120 x R120 block of the same pairs),
pred_C.npz (88 complex characters -> 176 one-sided L-functions chi, conj chi; 60 real by a seeded draw from R120;
all CC and CR pairs), onepoint.npz (one-point laws and resonance predictors), and freeze_summary.json.
usage: freeze.py STEP   (STEP = real | complex | onepoint | summary; each under 10 minutes)"""
import sys, math, os, json, time, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hlib import *

HERE = os.path.dirname(os.path.abspath(__file__))
ZD = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')) + '/'
P = 30000
SIG = (0.5, 1.0)
SEED_CR = 20261125
step = sys.argv[1]
t0 = time.time()


def keys(f):
    with np.load(ZD + f) as z:
        return list(z.files)


R120 = sorted(int(k[1:]) for f in ('r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz') for k in keys(f))
B = sorted(int(k[1:]) for k in keys('r134_B.npz'))
DN = sorted(int(k[1:]) for k in keys('r136_Dn.npz'))
C134 = sorted(int(k[1:]) for k in keys('r134_C.npz'))
OO = R120 + B + DN
CPLX = sorted((int(k.split('_')[0][1:]), int(k.split('_')[1][1:])) for i in range(4) for k in keys('r121_part%d.npz' % i))
assert len(OO) == 546 and len(set(OO)) == 546 and set(C134) == set(R120) and len(CPLX) == 88
rng = np.random.default_rng(SEED_CR)
CR_REAL = sorted(int(x) for x in rng.choice(np.array(R120), 60, replace=False))
pari = pari_instance()

FN = sorted(int(k[1:]) for k in keys('r138_Fn.npz'))       # USED set: validation design only
if step in ('real', 'fn'):
    LST, fname = (OO, 'pred_OO.npz') if step == 'real' else (FN, 'pred_FN.npz')
    pr, X = real_table(pari, LST, P)
    print('tables %.0fs' % (time.time() - t0), flush=True)
    Jm, H0 = J_real_all(X, pr, SIG, verbose=True)
    I, J = np.triu_indices(len(LST), 1)
    kap = np.array([phi_ratio(abs(D)) for D in LST])
    out = dict(Ds=np.array(LST), I=I, J=J, kappa=kap, H0=H0[I, J], P=P, umax=UMAX, du=DU)
    for s in SIG:
        Jp = Jm[s][I, J]
        out['J_%.1f' % s] = Jp
        out['pred_%.1f' % s] = kap[I] * kap[J] * Jp / (4 * math.pi ** 2)
        out['pred1_%.1f' % s] = kap[I] * kap[J] * s * math.sqrt(2 * math.pi) / (4 * math.pi ** 2) * np.ones(len(I))
        print('sigma %.1f: J range %.3f .. %.3f, mean %.4f; pred mean %.4f sd %.4f; symmetric check %.1e' % (
            s, Jp.min(), Jp.max(), Jp.mean(), out['pred_%.1f' % s].mean(), out['pred_%.1f' % s].std(), np.abs(Jm[s] - Jm[s].T).max()))
    np.savez_compressed(os.path.join(HERE, fname), **out)
    print('saved %s %.0fs' % (fname, time.time() - t0))

elif step == 'complex':
    # one-sided L-functions: index 2k = chi (zeros gamma > 0), 2k+1 = conj chi (zeros -gamma, gamma < 0)
    labs, qs, As, eps, eps_tau, tabs, kinds = [], [], [], [], [], [], []
    pr = primes_upto(P)
    for N, m in CPLX:
        q, a, e, et = conrey_info(pari, N, m)
        assert q == N, (N, m, q)
        _, c = conrey_table(pari, N, m, P)
        for cj in (0, 1):
            labs.append('N%d_m%d%s' % (N, m, '_conj' if cj else '')); qs.append(q); As.append(a)
            eps.append(np.conj(e) if cj else e); eps_tau.append(np.conj(et) if cj else et)
            tabs.append(np.conj(c) if cj else c); kinds.append(1 + cj)
    _, XR = real_table(pari, CR_REAL, P)
    for D, row in zip(CR_REAL, XR):
        labs.append('D%d' % D); qs.append(abs(D)); As.append(0 if D > 0 else 1); eps.append(1 + 0j); eps_tau.append(1 + 0j)
        tabs.append(row.astype(complex)); kinds.append(0)
    tabs = np.array(tabs); eps = np.array(eps); eps_tau = np.array(eps_tau); qs = np.array(qs); As = np.array(As); kinds = np.array(kinds)
    print('root numbers: max |eps_PARI - tau/(i^a sqrt q)| = %.1e; tables %.0fs' % (np.abs(eps - eps_tau).max(), time.time() - t0), flush=True)
    Jm = J_general_all(tabs, pr, SIG, P0=100, K=3, verbose=True)
    n = len(labs); nc = 176
    print('J done %.0fs; hermitian check max |J_ij - conj J_ji| = %.1e' % (time.time() - t0, max(np.abs(Jm[s] - Jm[s].conj().T).max() for s in SIG)), flush=True)
    I, J = np.triu_indices(nc, 1)                    # CC: all pairs among the 176 one-sided functions
    Ic = np.repeat(np.arange(nc), 60); Jr = np.tile(np.arange(nc, n), nc)   # CR
    kap = np.array([phi_ratio(int(q)) for q in qs])
    out = dict(labels=np.array(labs), q=qs, a=As, eps=eps, kappa=kap, kind=kinds, I_cc=I, J_cc=J, I_cr=Ic, J_cr=Jr, P=P,
               cr_real=np.array(CR_REAL), seed_cr=SEED_CR)
    for s in SIG:
        for nm, (a_, b_) in (('cc', (I, J)), ('cr', (Ic, Jr))):
            Jp = Jm[s][a_, b_]
            out['J_%s_%.1f' % (nm, s)] = Jp
            out['pred_%s_%.1f' % (nm, s)] = predict(kap[a_], kap[b_], eps[a_], eps[b_], Jp)
            out['pred1_%s_%.1f' % (nm, s)] = predict(kap[a_], kap[b_], eps[a_], eps[b_], s * math.sqrt(2 * math.pi))
            out['prednoroot_%s_%.1f' % (nm, s)] = kap[a_] * kap[b_] * Jp / (4 * math.pi ** 2)
            out['predconjJ_%s_%.1f' % (nm, s)] = predict(kap[a_], kap[b_], eps[a_], eps[b_], np.conj(Jp))
    # cross-check of the two code paths: real x real pairs among the 60 against J_real_all
    pr2, X60 = real_table(pari, CR_REAL, P)
    Jr60, _ = J_real_all(X60, pr2, SIG)
    for s in SIG:
        d = np.abs(Jm[s][nc:, nc:] - Jr60[s]); print('sigma %.1f: general vs real code on 60 x 60 real pairs: max |dJ| %.1e' % (s, d.max()))
    np.savez_compressed(os.path.join(HERE, 'pred_C.npz'), **out)
    print('saved pred_C.npz %.0fs' % (time.time() - t0))

elif step == 'onepoint':
    def resonance(chi_p, pr, q, a, w, nmax):
        tot = 0j
        for k in range(1, 40):
            m = pr.astype(float) ** k <= nmax
            if not m.any(): break
            n = pr[m].astype(np.int64) ** k
            tot += np.sum(chi_p[m] ** k * np.log(pr[m]) * np.exp(-2j * math.pi * (n % q) / q) * w(2 * math.pi * n / q))
        return -(1j ** a / math.sqrt(q)) * tot

    def windows(T1, T2):
        sharp = lambda t: ((np.asarray(t) > T1) & (np.asarray(t) <= T2)).astype(float)
        def hann(t):
            x = (np.asarray(t, float) - T1) / (T2 - T1); return np.where((x > 0) & (x <= 1), np.sin(math.pi * x) ** 2, 0.0)
        return sharp, hann

    out = {}
    for des, fl, T1, T2 in (('OO', [('R', D) for D in OO], 100.0, 600.0), ('HI', [('R', D) for D in R120], 600.0, 1000.0),
                             ('C', [('C', N, m, cj) for N, m in CPLX for cj in (0, 1)], 100.0, 500.0),
                             ('FN', [('R', D) for D in FN], 100.0, 600.0)):
        sharp, hann = windows(T1, T2)
        rows = []
        for f in fl:
            if f[0] == 'R':
                D = f[1]; q = abs(D); a = 0 if D > 0 else 1; e = 1 + 0j
                nmax = int(q * T2 / (2 * math.pi)) + 2
                pr = primes_upto(nmax)
                chi = np.array([int(x) for x in pari('apply(p->kronecker(%d,p), primes([2,%d]))' % (D, nmax))], float).astype(complex)
                lab = 'D%d' % D
            else:
                N, m, cj = f[1], f[2], f[3]
                q, a, e, _ = conrey_info(pari, N, m); nmax = int(q * T2 / (2 * math.pi)) + 2
                pr, chi = conrey_table(pari, N, m, nmax)
                if cj: chi = np.conj(chi); e = np.conj(e)
                lab = 'N%d_m%d%s' % (N, m, '_conj' if cj else '')
            kap = phi_ratio(q)
            M = -e * kap * (T2 - T1) / (2 * math.pi)
            rows.append((lab, q, a, kap, e, M, M / 2, resonance(chi, pr, q, a, sharp, nmax), resonance(chi, pr, q, a, hann, nmax)))
        out[des + '_labels'] = np.array([r[0] for r in rows]); out[des + '_q'] = np.array([r[1] for r in rows])
        out[des + '_a'] = np.array([r[2] for r in rows]); out[des + '_kappa'] = np.array([r[3] for r in rows])
        out[des + '_eps'] = np.array([r[4] for r in rows]); out[des + '_law_sharp'] = np.array([r[5] for r in rows])
        out[des + '_law_hann'] = np.array([r[6] for r in rows]); out[des + '_res_sharp'] = np.array([r[7] for r in rows])
        out[des + '_res_hann'] = np.array([r[8] for r in rows]); out[des + '_window'] = np.array([T1, T2])
        L = out[des + '_law_hann']; R = out[des + '_res_hann']
        print('%s: %d functions; resonance/law (Hann) complex slope %.4f%+.4fi; sd of (R - law)/|law| %.3f; %.0fs' % (
            des, len(rows), *(lambda c: (c.real, c.imag))(np.sum(np.conj(L) * R) / np.sum(np.abs(L) ** 2)), np.std(np.abs(R - L) / np.abs(L)), time.time() - t0), flush=True)
    np.savez_compressed(os.path.join(HERE, 'onepoint.npz'), **out)

elif step == 'summary':
    h = {}
    for f in ('pred_OO.npz', 'pred_C.npz', 'pred_FN.npz', 'onepoint.npz', 'hlib.py', 'freeze.py'):
        h[f] = hashlib.sha256(open(os.path.join(HERE, f), 'rb').read()).hexdigest()
    json.dump({'sha256': h, 'R120': R120, 'B': B, 'Dn': DN, 'complex': CPLX, 'cr_real': CR_REAL, 'seed_cr': SEED_CR, 'P': P,
               'sigmas': SIG, 'u_grid': [0, UMAX, DU]}, open(os.path.join(HERE, 'freeze_manifest.json'), 'w'), indent=1)
    for k, v in h.items(): print(v, k)
