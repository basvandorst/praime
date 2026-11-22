"""Round 59 (preregister/round59.txt): the distances d between consecutive twin-prime pairs, from the conditioned
expansion applied to twin events. y = 13: classes a mod 30030 with a, a+2 coprime to 30030 (1485). In class a an
interior twin start h (multiple of 6 in (0, d)) is possible iff a+h and a+h+2 are rough; it is a twin with chance
rho0^2 s_h; Pi_a (no interior twin) is expanded in Bahadur form to order 2. S_{>13} uses the exact residue counts for
every prime 17 <= q <= 3000 (all spans here are below 1210) and a generic tail beyond.
usage: r59.py law | data | compare"""
import os, sys, json, math, time
import numpy as np
import mpmath
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
PATH = os.path.join(OUT, 'r59.json' if int(os.environ.get('R59_Y', 13)) == 13 else 'r59_y%s.json' % os.environ.get('R59_Y'))
Y = int(os.environ.get('R59_Y', 13))
SMALL = [q for q in (2, 3, 5, 7, 11, 13, 17, 19) if q <= Y]; PP = math.prod(SMALL); DEN = math.prod(1 - 1 / q for q in SMALL)
MID = [int(q) for q in C.small_primes(3001) if q > Y]
BIGP = C.small_primes(20_000_000).astype(float); BIGP = BIGP[BIGP > 3000]
TAIL = {k: float(np.sum(np.log1p(-k / BIGP) - k * np.log1p(-1 / BIGP))) for k in (2, 4, 6, 8)}
L0, L1 = 10 ** 9, 2 * 10 ** 9
LAMS = [math.log(L0), math.log(1.5e9), math.log(L1)]
DMAX = 1200
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1 and math.gcd(a + 2, PP) == 1], dtype=np.int64)


def lnS_parts(d):
    """ln S_{>13} of B = {0,2,d,d+2} (minus its (1-1/q)^-4 normalisation handled in the ratios), ln s_h, and r_hk."""
    h = np.arange(6, d, 6, dtype=np.int64); n = len(h)
    lnB = 0.0; lns = np.zeros(n); lnr = np.zeros((n, n))
    for q in MID:
        B = {0 % q, 2 % q, d % q, (d + 2) % q}; nuB = len(B)
        Bm = np.zeros(q, bool); Bm[list(B)] = True
        r1 = h % q; r2 = (h + 2) % q; o1 = ~Bm[r1]; o2 = ~Bm[r2]
        c = o1.astype(int) + o2.astype(int)
        lnB += math.log(1 - nuB / q) - 4 * math.log(1 - 1 / q)
        nh = nuB + c
        lns += np.log1p(-nh / q) - math.log(1 - nuB / q) - 2 * math.log(1 - 1 / q)
        if n > 1:
            # coincidences between the new residues of k and of h
            delta = ((r1[None, :] == r1[:, None]) & o1[None, :]).astype(int) + ((r1[None, :] == r2[:, None]) & o1[None, :]) \
                  + ((r2[None, :] == r1[:, None]) & o2[None, :]) + ((r2[None, :] == r2[:, None]) & o2[None, :])
            nhk = nuB + c[:, None] + c[None, :] - delta
            with np.errstate(divide='ignore'):
                lnr += np.log1p(-nhk / q) + math.log(1 - nuB / q) - np.log1p(-nh / q)[:, None] - np.log1p(-nh / q)[None, :]
    lnB += TAIL[4]; lns += TAIL[6] - TAIL[4]; lnr += TAIL[8] + TAIL[4] - 2 * TAIL[6]
    np.fill_diagonal(lnr, 0.0)
    return h, lnB, lns, lnr


def share(d, lams):
    """P(next twin start at distance d | twin at 0), for each lam."""
    h, lnB, lns, lnr = lnS_parts(d)
    lnS02 = sum(math.log(1 - 2 / q) - 2 * math.log(1 - 1 / q) for q in MID) + TAIL[2]
    cls = CLS[(np.gcd(CLS + d, PP) == 1) & (np.gcd(CLS + d + 2, PP) == 1)]
    if len(cls) == 0:
        return [0.0] * len(lams)
    M = ((np.gcd(cls[:, None] + h[None, :], PP) == 1) & (np.gcd(cls[:, None] + h[None, :] + 2, PP) == 1)).astype(float) if len(h) else np.zeros((len(cls), 0))
    R = np.expm1(lnr)
    out = []
    for lam in lams:
        rho0 = 1 / (lam * DEN)
        p = rho0 ** 2 * np.exp(lns)
        L = M @ np.log1p(-p) if len(h) else np.zeros(len(cls))
        if len(h) > 1:
            w = p / (1 - p)
            L = L + 0.5 * np.sum((M @ (w[:, None] * w[None, :] * R)) * M, axis=1)
        out.append(float(np.sum(rho0 ** 2 * math.exp(lnB - lnS02) * np.exp(L)) / len(CLS)))
    return out


if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {'law': {}}
    mode = sys.argv[1]; t0 = time.time()
    if mode == 'law':
        for d in range(6, DMAX + 1, 6):
            if str(d) in res['law'] or time.time() - t0 > 420:
                continue
            res['law'][str(d)] = share(d, LAMS)
            json.dump(res, open(PATH, 'w'))
        print('law done', len(res['law']), 'of', DMAX // 6, round(time.time() - t0), 's', flush=True)
    elif mode == 'data':
        cnt = {}; last = None; n = 0
        for s in range(L0, L1, 10 ** 8):
            P = np.nonzero(C.primes_mask(s, 10 ** 8))[0].astype(np.int64) + s
            if last is not None:
                P = np.concatenate([[last], P])
            last = int(P[-1])
            tw = P[:-1][np.diff(P) == 2]
            cnt.setdefault('starts', []).extend(tw.tolist())
        st = np.array(cnt['starts']); dd = np.diff(st)
        v, c = np.unique(dd, return_counts=True)
        res['data'] = {str(a): int(b) for a, b in zip(v.tolist(), c.tolist())}; res['twins'] = len(st); res['gaps'] = len(dd)
        json.dump(res, open(PATH, 'w')); print('twins', len(st), 'gaps', len(dd), 'not multiple of 6:', int(np.sum(dd % 6 != 0)), flush=True)
    else:
        N = res['gaps']; data = {int(k): v for k, v in res['data'].items()}
        u = np.linspace(LAMS[0], LAMS[-1], 401); wt = np.exp(u) / u ** 2
        law = {}
        for k, v in res['law'].items():
            ls = np.interp(u, LAMS, v); law[int(k)] = float(np.sum(ls * wt) / np.sum(wt))
        mass = sum(law.values())
        lam_m = float(np.sum(u * wt) / np.sum(wt)); pi6 = 6 * 1.3203236316937 / lam_m ** 2
        bands = [(lo, lo + 24) for lo in range(6, 600, 30)] + [(lo, lo + 54) for lo in range(606, 1200, 60)]
        rows = []; chi = 0.0; dof = 0; chiP = 0.0
        for lo, hi in bands:
            c = sum(v for d, v in data.items() if lo <= d <= hi); e = N * sum(law.get(d, 0.0) for d in range(lo, hi + 1, 6))
            ep = N * sum(pi6 * (1 - pi6) ** (d // 6 - 1) for d in range(lo, hi + 1, 6))
            rows.append({'band': [lo, hi], 'count': c, 'law': e, 'poisson': ep, 'ratio': c / e if e else None, 'z': (c - e) / math.sqrt(e) if e else None, 'z_poisson': (c - ep) / math.sqrt(ep)})
            if e >= 50:
                chi += (c - e) ** 2 / e; chiP += (c - ep) ** 2 / ep; dof += 1
        tail_c = sum(v for d, v in data.items() if d > 1200); tail_e = N * (1 - mass)
        res['compare'] = {'rows': rows, 'chi2': chi, 'dof': dof, 'chi2_poisson': chiP, 'law_mass_to_1200': mass,
                          'tail_count': tail_c, 'tail_law': tail_e, 'mean_lambda': lam_m}
        for r in rows:
            print(r['band'], r['count'], round(r['law'], 1), round(r['poisson'], 1), r['ratio'] and round(r['ratio'], 4), r['z'] and round(r['z'], 2), round(r['z_poisson'], 1))
        print('chi2 law %.1f / %d dof; Poisson %.1f; tail >1200: %d vs %.1f; mass %.6f' % (chi, dof, chiP, tail_c, tail_e, mass), flush=True)
        json.dump(res, open(PATH, 'w'), indent=1)
