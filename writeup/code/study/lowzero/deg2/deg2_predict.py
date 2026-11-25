"""Frozen per-curve predictions of Sigma2(l) for a curve set (deg2lib.py), under the law and the alternatives, averaged
over the band with the zero-density weight (10 Gauss nodes in t). One cache file per curve (resumable); 'pack' writes the
predictions JSON (l = 0.5, 1, ..., 10) and prints its sha256.
  full    the law (ratios recipe, all local factors)
  k1      traversals k >= 2 of every good prime carry the Sato-Tate mean weight 2, not (alpha^k + beta^k)^2
  offzeta off-diagonal amplitude of degree-1 form |zeta(1+ie)|^2 Phi_off prod_{p|N} F_p (no Sym2, no a_p); diagonal = law
  deg1    every good prime has weight 1 at every traversal, bad primes missing (BK 2013 with the conductor's primes)
  gue     Sigma2 of GUE
  dP_K    the diagonal Sigma2 term of traversal K (1, 2, >= 3) of p = P in (2, 3, 5, 7) with its weight W_p(k)
  k1pP    the law with the traversals k >= 2 of P (2 or 3) at the Sato-Tate mean weight 2 (template of the doubled orbit)
  w1pP    the law with the first traversal of P at weight 1 (template of the first-traversal weight a_p^2/p)
usage: deg2_predict.py SET.json CACHEDIR PART NPARTS [SECONDS]  |  deg2_predict.py SET.json CACHEDIR pack OUT.json"""
import sys, os, json, math, time, pickle, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deg2lib as D

LS = [0.5 * i for i in range(1, 21)]
IDX = [int(round(l / D.W)) - 1 for l in LS]
SPLIT = (2, 3, 5, 7)


def weights(c, p):
    """W_p(k), k = 1..K for the curve c (good: (2 cos k theta)^2; multiplicative: p^-k)"""
    K = 60
    if p in c.bad:
        return [p ** (-k) for k in range(1, K + 1)]
    lam = c.ap[p] / math.sqrt(p); th = math.acos(max(-1.0, min(1.0, lam / 2)))
    return [(2 * math.cos(k * th)) ** 2 for k in range(1, K + 1)]


def diag_pieces(N, bad, ap, p, lo, hi, mode='law'):
    """height-averaged diagonal Sigma2 terms of prime p, traversals k = 1, 2, >= 3 (closed form).
    mode 'law': W_p(k); 'st': W_p(1) as the law, W_p(k >= 2) = 2 (Sato-Tate mean); 'one': W_p(1) = 1, k >= 2 as the law."""
    Ls, ws = D.L_nodes(N, lo, hi, 10)
    if p in bad: Wk = [p ** (-k) for k in range(1, 61)]
    else:
        lam = ap / math.sqrt(p); th = math.acos(max(-1.0, min(1.0, lam / 2)))
        Wk = [(2 * math.cos(k * th)) ** 2 for k in range(1, 61)]
        if mode == 'st': Wk = [Wk[0]] + [2.0] * 59
        if mode == 'one': Wk = [1.0] + Wk[1:]
    out = {}
    for kk, ks in (('1', [1]), ('2', [2]), ('3', range(3, 61))):
        out[kk] = sum(w * sum(D.dsig_diag_k(p, Wk[k - 1], k, L, D.LGRID) for k in ks if p ** (-k) > 1e-18) for L, w in zip(Ls, ws))
    return out


def nz_expected(N, lo, hi):
    return float((D.theta_E(N, hi) - D.theta_E(N, lo)) / math.pi)


def predict(Z, c, lo, hi):
    Ls, ws = D.L_nodes(c.N, lo, hi, 10)
    assert 2 * math.pi * D.LMAX / min(Ls) <= D.EMAX, 'eps grid too short'
    g = D.sigma2_gue(); out = {}
    for mode in ('full', 'k1', 'offzeta'):
        out[mode] = sum(w * (g + D.cum_sigma(D.C_full(Z, c, L, kmode=mode))) for L, w in zip(Ls, ws))
    out['deg1'] = sum(w * (g + D.cum_sigma(D.C_deg1(Z, c.bad, L))) for L, w in zip(Ls, ws))
    out['gue'] = g
    for p in SPLIT:
        Wk = weights(c, p)
        for kk, ks in (('1', [1]), ('2', [2]), ('3', range(3, 61))):
            out['d%d_%s' % (p, kk)] = sum(w * sum(D.dsig_diag_k(p, Wk[k - 1], k, L, D.LGRID) for k in ks if p ** (-k) > 1e-18)
                                         for L, w in zip(Ls, ws))
    out['meanL'] = float(np.dot(ws, Ls))
    return out


def main():
    S = json.load(open(sys.argv[1])); cache = sys.argv[2]; os.makedirs(cache, exist_ok=True)
    lo, hi = S['band']
    if sys.argv[3] == 'pack':
        rows = []
        for c in S['curves']:
            f = os.path.join(cache, c['id'] + '.pkl')
            o = pickle.load(open(f, 'rb'))
            row = {'id': c['id'], 'ainv': c['ainv'], 'N': c['N'], 'w': c['w'], 'g2': c['g2'], 'g3': c['g3'],
                   'ap': c['ap'], 'bad': c['bad'], 'meanL': o['meanL'], 'LS1': o['LS1']}
            for p in (2, 3):                                     # per-prime templates (closed-form diagonal terms)
                law = diag_pieces(c['N'], c['bad'], c['ap'][str(p)], p, lo, hi, 'law')
                st = diag_pieces(c['N'], c['bad'], c['ap'][str(p)], p, lo, hi, 'st')
                one = diag_pieces(c['N'], c['bad'], c['ap'][str(p)], p, lo, hi, 'one')
                assert np.allclose(law['1'] + law['2'] + law['3'], o['d%d_1' % p] + o['d%d_2' % p] + o['d%d_3' % p], atol=1e-12)
                o['k1p%d' % p] = o['full'] - (law['2'] + law['3']) + (st['2'] + st['3'])   # repetitions of p at weight 2
                o['w1p%d' % p] = o['full'] - law['1'] + one['1']                          # first traversal of p at weight 1
            for k, v in o.items():
                if isinstance(v, np.ndarray): row[k] = [round(float(v[i]), 7) for i in IDX]
            row['nz_expected'] = round(nz_expected(c['N'], lo, hi), 2)
            rows.append(row)
        J = {'band': [lo, hi], 'l': LS, 'law': 'deg2lib.py (ratios recipe for L(s,E); see derivation.md)',
             'gue': [round(float(D.sigma2_gue()[i]), 7) for i in IDX], 'curves': rows}
        txt = json.dumps(J, sort_keys=True, separators=(',', ':'))
        open(sys.argv[4], 'w').write(txt)
        print('%d curves packed; sha256 %s' % (len(rows), hashlib.sha256(txt.encode()).hexdigest()))
        return
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
    zf = os.path.join(cache, 'zgrid.pkl')
    Z = pickle.load(open(zf, 'rb')) if os.path.exists(zf) else D.ZetaGrid()
    if not os.path.exists(zf): pickle.dump(Z, open(zf, 'wb'))
    part, nparts = int(sys.argv[3]), int(sys.argv[4]); budget = float(sys.argv[5]) if len(sys.argv) > 5 else 540
    t0 = time.time(); done = 0
    for c in S['curves'][part::nparts]:
        f = os.path.join(cache, c['id'] + '.pkl')
        if os.path.exists(f): continue
        if time.time() - t0 > budget: break
        t1 = time.time(); cv = D.Curve(pari, c['ainv']); assert cv.N == c['N']
        o = predict(Z, cv, lo, hi); o['LS1'] = cv.LS1
        pickle.dump(o, open(f, 'wb')); done += 1
        print('%s N=%d L=%.2f full(1,4,8)=%.4f %.4f %.4f  %.0fs' % (c['id'], c['N'], o['meanL'], o['full'][IDX[1]], o['full'][IDX[7]],
                                                                o['full'][IDX[15]], time.time() - t1), flush=True)
    left = sum(not os.path.exists(os.path.join(cache, c['id'] + '.pkl')) for c in S['curves'][part::nparts])
    print('part %d/%d: computed %d, left %d, %.0fs' % (part, nparts, done, left, time.time() - t0))


if __name__ == '__main__':
    main()
