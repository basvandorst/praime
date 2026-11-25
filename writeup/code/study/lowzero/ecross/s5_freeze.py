"""Step 5: per-pair predictions of the recipe (H1) and its alternatives, frozen with sha256.
   pred      = w1 w2 c1 c2 J/(4 pi^2)            (H1: the recipe's double-swap term)
   pred1     = w1 w2 c1 c2 sigma sqrt(2 pi)/(4 pi^2)   (H = 1 part: the product of the two one-point laws)
   prednoroot= c1 c2 J/(4 pi^2)                  (H1 without root numbers)
   predrev   = w1 w2 c1 c2 Jrev/(4 pi^2)         (lines paired r with 1/r: reversed direction)
   H0: 0;  H1 with the wrong sign: -pred.
J = int H(u) K(u) du, K = exp(-u^2/2 sigma^2), sigma = 0.5 and 1.0; H from exact local factors p <= 3000 and first-order
lines (+ the non-oscillating second-order term) 3000 < p <= 1e5; u grid 0..6 step 0.05 (H(-u) = conj H(u)).
usage: s5_freeze.py check | freeze"""
import sys, json, math, hashlib, time
import numpy as np
sys.path.insert(0, '.')
import elib as L
C = np.load('curves.npz')
ids = [str(s) for s in C['ids']]; sets = [str(s) for s in C['set']]
U = np.arange(0, 6.0001, 0.05); SIGMAS = (0.5, 1.0)
def cur(i):
    return dict(primes=C['primes'], ap=C['ap'][i].astype(np.int64), bad=json.loads(str(C['bad'][i])))
mode = sys.argv[1]
LINES = [(1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (4, 1), (1, 4), (5, 1), (1, 5), (7, 1), (1, 7), (3, 2), (2, 3)]
def rho_r(i, a, b):
    bad = json.loads(str(C['bad'][i])); out = 1.0
    for n_, sgn in ((a, 1), (b, -1)):
        m = n_; q = 2
        while m > 1:
            if m % q == 0:
                v = 0
                while m % q == 0: m //= q; v += 1
                k = int(np.searchsorted(C['primes'], q))
                out *= float(L.rho_pv(q, np.array([float(C['ap'][i][k])]), np.array([q in bad]), sgn * v)[0])
            q += 1
    return out
if mode == 'check':
    out = open('s5_check_out.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
    sel = [i for i in range(len(ids)) if sets[i] == 'dev']
    cs = [cur(i) for i in sel]
    t0 = time.time()
    ref = L.lnH_matrix(cs, U, 3000, 10 ** 5); P('dev 20x20, P0 3000, P 1e5: %.0fs' % (time.time() - t0))
    Jref = {s: L.J_from_H(np.exp(ref), U, s) for s in SIGMAS}
    iu = np.triu_indices(len(sel), 1)
    for s in SIGMAS:
        J = Jref[s][iu]; P('sigma %.1f: J/(sigma sqrt(2pi)) over the 190 dev pairs: min %.3f median %.3f max %.3f' % (s, (J / (s * math.sqrt(2 * math.pi))).min(), np.median(J / (s * math.sqrt(2 * math.pi))), (J / (s * math.sqrt(2 * math.pi))).max()))
    for (P0, PP, sec) in ((3000, 10 ** 4, True), (3000, 3 * 10 ** 4, True), (1000, 10 ** 5, True), (300, 10 ** 5, True), (100, 10 ** 5, True), (3000, 10 ** 5, False)):
        t = L.lnH_matrix(cs, U, P0, PP, second=sec)
        for s in SIGMAS:
            d = (L.J_from_H(np.exp(t), U, s) - Jref[s])[iu]
            cc = (C['c'][sel][:, None] * C['c'][sel][None, :])[iu] / (4 * math.pi ** 2)
            P('  P0 %5d P %6d second %d, sigma %.1f: max |dJ| %.2e, max |dJ|/|J| %.2e, max |d pred| %.2e' % (P0, PP, sec, s, np.abs(d).max(), (np.abs(d) / np.abs(Jref[s][iu])).max(), np.abs(d * cc).max()))
    # quadrature: step 0.025 to 8
    U2 = np.arange(0, 8.0001, 0.025); t = L.lnH_matrix(cs, U2, 3000, 10 ** 5)
    for s in SIGMAS:
        d = (L.J_from_H(np.exp(t), U2, s) - Jref[s])[iu]; P('  u-grid step 0.025 to 8, sigma %.1f: max |dJ| %.2e' % (s, np.abs(d).max()))
    for s in SIGMAS:
        cc = (C['c'][sel][:, None] * C['c'][sel][None, :])[iu] / (4 * math.pi ** 2)
        P('sigma %.1f: dev pred = w w c c J/4pi^2: min %.4f median %.4f max %.4f; c1 c2 J/(4 pi^2 sigma sqrt 2pi): min %.3f max %.3f' % (s, (cc * Jref[s][iu] * (C['w'][sel][:, None] * C['w'][sel][None, :])[iu]).min(), np.median(cc * Jref[s][iu] * (C['w'][sel][:, None] * C['w'][sel][None, :])[iu]), (cc * Jref[s][iu] * (C['w'][sel][:, None] * C['w'][sel][None, :])[iu]).max(), (cc * Jref[s][iu]).min() / (s * math.sqrt(2 * math.pi)), (cc * Jref[s][iu]).max() / (s * math.sqrt(2 * math.pi))))
    # lines form against the product: J = sigma sqrt(2pi) sum_r h_r exp(-sigma^2 ln^2 r/2), r over p^v, p <= 13 only, for a few pairs
    small = np.log(np.prod([L.h_lines(p, np.array([float(c_['ap'][int(np.searchsorted(C['primes'], p))]) for c_ in cs[:4]]), np.array([p in c_['bad'] for c_ in cs[:4]]), np.array([float(c_['ap'][int(np.searchsorted(C['primes'], p))]) for c_ in cs[:4]]), np.array([p in c_['bad'] for c_ in cs[:4]]), 1j * U) for p in (2, 3, 5)], axis=0))
    for s in SIGMAS:
        Jp = L.J_from_H(np.exp(small), U, s)
        # direct: sum over v2,v3,v5,v7,v11,v13 of prod rho rho p^-|v| exp(-s^2 (sum v ln p)^2/2)
        import itertools
        prs = [2, 3, 5]; rng = {2: range(-2, 34), 3: range(-2, 24), 5: range(-2, 18)}
        mx = 0
        for a_ in range(4):
            for b_ in range(4):
                fac = []
                for p in prs:
                    k = int(np.searchsorted(C['primes'], p))
                    r1 = {v: float(L.rho_pv(p, np.array([float(cs[a_]['ap'][k])]), np.array([p in cs[a_]['bad']]), v)[0]) for v in rng[p]}
                    r2 = {v: float(L.rho_pv(p, np.array([float(cs[b_]['ap'][k])]), np.array([p in cs[b_]['bad']]), v)[0]) for v in rng[p]}
                    fac.append([(v * math.log(p), r1[v] * r2[v] * p ** (-abs(v))) for v in rng[p] if abs(r1[v] * r2[v] * p ** (-abs(v))) > 1e-16])
                tot = 0.0
                for combo in itertools.product(*fac):
                    lr = sum(c_[0] for c_ in combo); hh = np.prod([c_[1] for c_ in combo])
                    if abs(hh) < 1e-14: continue
                    tot += hh * math.exp(-s * s * lr * lr / 2)
                tot *= s * math.sqrt(2 * math.pi)
                mx = max(mx, abs(tot - Jp[a_, b_]))
        P('  lines sum (r = 2^a 3^b 5^c) against int h_2 h_3 h_5 K (quadrature), 16 pairs, sigma %.1f: max diff %.1e' % (s, mx))
    out.close()
else:
    out = open('s5_freeze_out.txt', 'w')
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
    man = {}
    for setname in ('dev', 'fresh'):
        sel = [i for i in range(len(ids)) if sets[i] == setname]
        cs = [cur(i) for i in sel]; t0 = time.time()
        lnH = L.lnH_matrix(cs, U, 3000, 10 ** 5); lnHr = L.lnH_matrix(cs, U, 3000, 10 ** 5, reverse=True)
        H = np.exp(lnH); Hr = np.exp(lnHr); del lnH, lnHr
        I, J = np.triu_indices(len(sel), 1)
        w = C['w'][sel].astype(float); c = C['c'][sel]
        D = dict(ids=np.array([ids[i] for i in sel]), N=C['N'][sel], w=C['w'][sel], c=c, LS1=C['LS1'][sel], B=C['B'][sel], I=I, J=J,
                 u=U, H0=H[I, J, 0].real)
        for s in SIGMAS:
            Js = L.J_from_H(H, U, s); Jr = L.J_from_H(Hr, U, s)
            base = w[I] * w[J] * c[I] * c[J] / (4 * math.pi ** 2)
            D['J_%.1f' % s] = Js[I, J]; D['Jrev_%.1f' % s] = Jr[I, J]
            D['pred_%.1f' % s] = base * Js[I, J]
            D['pred1_%.1f' % s] = base * s * math.sqrt(2 * math.pi)
            D['prednoroot_%.1f' % s] = c[I] * c[J] / (4 * math.pi ** 2) * Js[I, J]
            D['predrev_%.1f' % s] = base * Jr[I, J]
            asym = np.abs(Js - Js.T).max()
            P('%s sigma %.1f: %d pairs; J/(sigma sqrt 2pi): min %.3f median %.3f max %.3f; pred: mean %.4f sd %.4f min %.4f max %.4f; pred1 sd %.4f; '
              'corr(pred, pred1) %.3f; corr(pred, predrev) %.3f; |J_ij - J_ji| max %.1e' % (
                setname, s, len(I), (Js[I, J] / (s * math.sqrt(2 * math.pi))).min(), np.median(Js[I, J] / (s * math.sqrt(2 * math.pi))),
                (Js[I, J] / (s * math.sqrt(2 * math.pi))).max(), D['pred_%.1f' % s].mean(), D['pred_%.1f' % s].std(),
                D['pred_%.1f' % s].min(), D['pred_%.1f' % s].max(), D['pred1_%.1f' % s].std(),
                np.corrcoef(D['pred_%.1f' % s], D['pred1_%.1f' % s])[0, 1], np.corrcoef(D['pred_%.1f' % s], D['predrev_%.1f' % s])[0, 1], asym))
        # one-point lines (per unit T/2pi): -w c rho(a/b)/sqrt(ab) for the registered list; synthetic-check data (p = 2, 3)
        D['lines'] = np.array(LINES)
        D['lineamp'] = np.array([[-C['w'][i] * C['c'][i] * rho_r(i, a, b) / math.sqrt(a * b) for (a, b) in LINES] for i in sel])
        D['lineamp_noroot'] = np.abs(C['w'][sel])[:, None] * 0 + np.array([[-C['c'][i] * rho_r(i, a, b) / math.sqrt(a * b) for (a, b) in LINES] for i in sel])
        for p, V in ((2, 40), (3, 26)):
            k = int(np.searchsorted(C['primes'], p)); vv = np.arange(-2, V + 1)
            D['M%d_v' % p] = vv
            D['M%d' % p] = np.array([[float(L.rho_pv(p, np.array([float(C['ap'][i][k])]), np.array([p in json.loads(str(C['bad'][i]))]), int(v))[0]) * p ** (-abs(v) / 2) for v in vv] for i in sel])
        H23 = np.ones((len(sel), len(sel), len(U)), complex)
        for p in (2, 3):
            k = int(np.searchsorted(C['primes'], p))
            a = C['ap'][sel][:, k].astype(float); b = np.array([p in json.loads(str(C['bad'][i])) for i in sel])
            H23 *= L.h_lines(p, a, b, a, b, 1j * U)
        for s in SIGMAS:
            D['J23_%.1f' % s] = L.J_from_H(H23, U, s)[I, J]
        f = 'pred_%s.npz' % setname
        np.savez(f, **D); man[f] = hashlib.sha256(open(f, 'rb').read()).hexdigest()
        P('%s frozen in %.0fs, sha256 %s' % (f, time.time() - t0, man[f]))
    for f in ('elib.py', 's5_freeze.py', 'curves.npz'):
        man[f] = hashlib.sha256(open(f, 'rb').read()).hexdigest()
    json.dump(man, open('freeze_manifest.json', 'w'), indent=1)
    out.close()
