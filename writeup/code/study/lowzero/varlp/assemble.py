"""assemble.py: the per-bad-set change of Var log|L'(rho)| at fixed L from the pieces (derivation.md sect. 4):
  dVar(S, L) = sum_{p in S} [ -Li2(1/p)/2 + dB_off + dA_diag + dA_off + dV3_Q + dV3_P ]
             + sum_{pairs in S} [ NA_B + NA_A + NA_V3P ]                 (non-additive parts; pairs among 2..13 only)
             - [ (mu + dmu(S))^2 - mu^2 ],   dmu(S) = sum_p (dmu_diag + dmu_off) + sum_pairs NA_mu,
with mu = mu(L) of the base (mu.py), the one-dimensional pieces from oned2.py, dV3_P and NA_V3P from twod2.py
(runs/t2_L*_V*.json). Pieces are interpolated in 1/L (cubic) across the computed L grid; primes not computed are
interpolated in ln p of p*piece. Importable: from assemble import Theory.
"""
import sys, os, json, glob, math
import numpy as np
from scipy.interpolate import CubicSpline
from mpmath import polylog
HERE = os.path.dirname(os.path.abspath(__file__))


class Theory:
    def __init__(self, runs=os.path.join(HERE, 'runs', 't2_L*_V70.json')):
        mu = json.load(open(os.path.join(HERE, 'mu.json')))
        Lm = np.array([r['L'] for r in mu['rows']]); self.mu_sp = CubicSpline(1 / Lm[::-1], np.array([r['mu'] for r in mu['rows']])[::-1])
        o = json.load(open(os.path.join(HERE, 'oned2.json'))); Lo = np.array(o['L'])
        self.P = sorted(int(p) for p in o['singles'])
        self.sp1 = {}
        for p in self.P:
            rows = o['singles'][str(p)]
            for k in ('dmu_diag', 'dmu_off', 'dA_diag', 'dA_off', 'dB_off', 'dV3_Q'):
                self.sp1[(p, k)] = CubicSpline(1 / Lo[::-1], np.array([r[k] for r in rows])[::-1])
        self.na1 = {}
        for key, rows in o['pairs'].items():
            pq = tuple(int(x) for x in key.split(','))
            for k in ('na_mu', 'na_A', 'na_B'): self.na1[(pq, k)] = CubicSpline(1 / Lo[::-1], np.array([r[k] for r in rows])[::-1])
        T = [json.load(open(f)) for f in sorted(glob.glob(runs))]; T.sort(key=lambda t: t['L'])
        self.L3 = np.array([t['L'] for t in T])
        self.sp3 = {}
        for p in self.P:
            y = np.array([[s['dV3P'] for s in t['singles'] if s['p'] == p][0] for t in T])
            self.sp3[p] = CubicSpline(1 / self.L3[::-1], y[::-1]) if len(T) > 2 else None
            sub = [[s for s in t['singles'] if s['p'] == p][0] for t in T]
            u = np.array([0.5 * (max(s['dV3P_sub'] + [s['dV3P']]) - min(s['dV3P_sub'] + [s['dV3P']])) for s in sub])
            self.sp3u = getattr(self, 'sp3u', {}); self.sp3u[p] = CubicSpline(1 / self.L3[::-1], np.abs(u[::-1]))
        self.na3 = {}
        for j, pr in enumerate(T[0]['pairs']):
            # pair non-additivity: mean over the four damping widths (its far field holds near-resonant terms
            # n ~ e^L that are below the conjecture's accuracy; their spread is carried as theory uncertainty)
            y = np.array([np.mean(t['pairs'][j]['na']) for t in T]); self.na3[tuple(pr['S'])] = CubicSpline(1 / self.L3[::-1], y[::-1])
            s = np.array([0.5 * (max(t['pairs'][j]['na'] + [t['pairs'][j]['na_ext']]) - min(t['pairs'][j]['na'] + [t['pairs'][j]['na_ext']])) for t in T])
            self.na3_unc = getattr(self, 'na3_unc', {}); self.na3_unc[tuple(pr['S'])] = CubicSpline(1 / self.L3[::-1], s[::-1])
        self.li2 = {}

    def _li2(self, p):
        if p not in self.li2: self.li2[p] = float(polylog(2, 1.0 / p))
        return self.li2[p]

    def _single_known(self, p, L):
        x = 1 / L
        d = {k: float(self.sp1[(p, k)](x)) for k in ('dmu_diag', 'dmu_off', 'dA_diag', 'dA_off', 'dB_off', 'dV3_Q')}
        d['dV3_P'] = float(self.sp3[p](x)); return d

    def single(self, p, L):
        """per-prime pieces at L (dict); p not on the computed list: p*piece interpolated linearly in ln p"""
        if p in self.P: d = self._single_known(p, L)
        else:
            big = [q for q in self.P if q >= 13]
            if p > big[-1]: raise ValueError('prime beyond the computed list: %d' % p)
            j = np.searchsorted(big, p); a, b = big[j - 1], big[j]
            da, db = self._single_known(a, L), self._single_known(b, L); t = (math.log(p) - math.log(a)) / (math.log(b) - math.log(a))
            d = {k: ((1 - t) * a * da[k] + t * b * db[k]) / p for k in da}
            lp = math.log(p); l1 = math.log(p / (p - 1)); d['dmu_diag'] = lp * l1 / L; d['dV3_Q'] = lp * l1 ** 2 / L
        d['dB_diag'] = -0.5 * self._li2(p)
        d['dmu'] = d['dmu_diag'] + d['dmu_off']
        d['lin'] = d['dB_diag'] + d['dB_off'] + d['dA_diag'] + d['dA_off'] + d['dV3_Q'] + d['dV3_P']
        return d

    def unc(self, S, L, floor=3e-4):
        """numerical uncertainty of dvar (not the conjecture's): damping spread of dV3_P per prime and of the pair terms,
        plus a floor for the domain size V (V = 60 vs 70: 1.4e-4 at L = 10) and the L interpolation; in quadrature."""
        S = sorted(set(S)); x = 1 / L; v = floor ** 2
        for p in S:
            q = p if p in self.P else min((r for r in self.P if r >= p), default=self.P[-1])
            v += (float(self.sp3u[q](x)) * (q / p if q != p else 1)) ** 2
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                if (S[i], S[j]) in self.na3_unc: v += float(self.na3_unc[(S[i], S[j])](x)) ** 2
        return math.sqrt(v)

    def dvar(self, S, L, parts=False):
        S = sorted(set(S)); mu = float(self.mu_sp(1 / L))
        lin = 0.0; dmu = 0.0; na = 0.0
        for p in S:
            d = self.single(p, L); lin += d['lin']; dmu += d['dmu']
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                pq = (S[i], S[j])
                if (pq, 'na_mu') in self.na1:
                    x = 1 / L
                    na += float(self.na1[(pq, 'na_A')](x) + self.na1[(pq, 'na_B')](x) + self.na3[pq](x)); dmu += float(self.na1[(pq, 'na_mu')](x))
        tot = lin + na - ((mu + dmu) ** 2 - mu ** 2)
        if parts: return dict(total=tot, lin=lin, na=na, dmu=dmu, mu=mu, musq=-((mu + dmu) ** 2 - mu ** 2))
        return tot


if __name__ == '__main__':
    th = Theory()
    print('computed L grid (three-point runs):', list(th.L3))
    print('\nper-prime change of Var log|L\'| at fixed L (bad prime added to S = {}):')
    print('  L     p   -Li2/2    dB_off   dA(diag+off)  dV3_Q    dV3_P   -2mu dmu - dmu^2  | total    first-order (LDM)')
    for L in (8.0, 9.0, 10.0, 10.2, 11.0, 11.4, 12.0):
        for p in (2, 3, 5, 7, 11, 13):
            d = th.single(p, L); f = th.dvar([p], L, parts=True); l1 = math.log(p / (p - 1))
            print('%5.1f %3d  %+.5f  %+.5f  %+.5f     %+.5f  %+.5f  %+.5f         | %+.5f  %+.5f' % (
                L, p, d['dB_diag'], d['dB_off'], d['dA_diag'] + d['dA_off'], d['dV3_Q'], d['dV3_P'], f['musq'], f['total'],
                d['dB_diag'] + math.log(p) * l1 * (2 + l1) / L))


def fit_orders(th, primes=(2, 3, 5, 7, 11, 13), Ls=np.arange(8.5, 12.01, 0.25)):
    """dVar(p, L) = -Li2(1/p)/2 + a1/L + a2/L^2 + a3/L^3: a1 exact (4.4); a2, a3 least squares over 8.5 <= L <= 12"""
    out = {}
    for p in primes:
        l1 = math.log(p / (p - 1)); a1 = math.log(p) * l1 * (2 + l1)
        y = np.array([th.dvar([p], L) + 0.5 * th._li2(p) - a1 / L for L in Ls])
        A = np.vstack([Ls ** -2.0, Ls ** -3.0]).T; c, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        out[p] = (a1, float(c[0]), float(c[1]), float(np.abs(A @ c - y).max()))
    return out


if __name__ == '__main__':
    print('\nexpansion dVar = -Li2(1/p)/2 + a1/L + a2/L^2 + a3/L^3 (a1 exact, a2 a3 fitted on 8.5 <= L <= 12; max misfit):')
    for p, (a1, a2, a3, mf) in fit_orders(th).items():
        print('  p=%2d  a1 = %.4f  a2 = %+.3f  a3 = %+.2f   (misfit %.1e)' % (p, a1, a2, a3, mf))
