#!/usr/bin/env python3
"""ecross: registered test of the swapped (Gram-lattice) cross term between the zeros of two elliptic-curve L-functions.

Statistic (per pair of curves i < j; the weights cancel the recipe's phase exactly):
   X_ij(sigma) = (1/T) sum_{g_i, g_j} e^{2i theta_i(g_i) - 2i theta_j(g_j)} K(g_j - g_i) w((g_i + g_j)/2),
   computed as (1/(sqrt(pi) s T)) int_W G_i(t) conj G_j(t) dt with G_i(t) = sum_g e^{2i theta_i(g)} exp(-(t - g)^2/2 s^2),
   s = sigma/sqrt 2, K(u) = exp(-u^2/2 sigma^2), W = [20, 396], T = 376, t-grid step 0.02, sigma = 0.5 (primary), 1.0.
   theta_E(t) = (t/2) ln N - t ln 2 pi + Im lnGamma(1 + i t) (scipy loggamma).
Law (H1, frozen in pred_<design>.npz):  E X_ij = (1/4 pi^2) w_i w_j c_i c_j J_ij, J_ij = int H_ij(u) K(u) du,
   H_ij(u) = prod_p sum_v rho_i(p^v) rho_j(p^v) p^{-|v|} p^{iuv} = |L(1+iu, E_i x E_j)|^2 A_ij(iu); c = B/L(1,Sym^2 E).
One-point sums: G_E(xi) = sum_{20 < g <= 400} e^{2i theta(g) - i xi g} at xi = ln r, law -w c rho(r)/sqrt(ab) (T2-T1)/2pi
   (frozen 'lineamp' times 380/2pi); r = 1 was computed on all curves before registration (not blind); r != 1 is blind
   on the fresh curves.
usage:  ecross_test.py STEP [--design fresh|dev] [--zeros PATH] [--frozen DIR] [--out DIR] [--synthetic poisson|gram|lines]
        [--seed S]       STEP = onepoint | x | readings | all   (cached in --out; resumable)
Completeness of the zero lists is assumed from round 144 (all 138 fresh sets complete after the deduplication; the 20
development sets complete).  Needs numpy, scipy only."""
import sys, os, math, json, time, hashlib, argparse
import numpy as np
from scipy.special import loggamma

FROZEN_SHA = {'pred_fresh.npz': 'ef45af2e11c905acbfb287dd8fd786322f1b29afe251d77e81f97369af765010',
              'pred_dev.npz': '25833cd6c2e9334db0fb9cc98e2a8235464dc4f42005ee4e3bf298bc4030c3e2'}      # filled at freezing; the test refuses other files
SIGMAS = (0.5, 1.0)
DT = 0.02
BAND = (20.0, 400.0)          # one-point band (T1, T2]
WIN = (20.0, 396.0)           # X window (zeros stored to >= 401: 6 s = 4.3 inside the top for sigma = 1)
CONTROL_XI = (0.35, 0.55, 0.9, 1.3, -0.35, -0.55, -0.9, -1.3)
DEFAULT_ZEROS = {'fresh': os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r144', 'fresh_zeros.npz'),
                 'dev': os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r148', 'dev_zeros.npz')}


def theta_E(N, t):
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(2 * math.pi) + np.imag(loggamma(1 + 1j * t))


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def load_frozen(fdir, design):
    f = 'pred_%s.npz' % design; p = os.path.join(fdir, f)
    if FROZEN_SHA.get(f) and sha(p) != FROZEN_SHA[f]:
        raise SystemExit('frozen file %s does not match its registered sha256' % f)
    return dict(np.load(p, allow_pickle=False))


def load_zeros(design, F, path):
    ids = [str(s) for s in F['ids']]
    if path.endswith('.npz'):
        with np.load(path) as z:
            return [np.sort(np.array(z[i], float)) for i in ids]
    return [np.sort(np.load(os.path.join(path, i + '.npy'))) for i in ids]


# ------------------------------------------------------------------ synthetic zero sets (dry runs)
def synth_zeros(F, kind, seed, top=402.0):
    """'poisson': independent points, density theta'/pi (no term); 'gram': one point per Gram interval at theta/pi =
    k + phi + j xi (phi = 1/2 for w = +1, 0 for w = -1; e^{-2 pi^2 j^2} = c/L(t)): the r = 1 line exactly (H = 1 part);
    'lines': Poisson with intensity theta'/pi + (1/pi) Re[a(t) e^{-2i theta(t)}], a(t) = -w c M_2(t) M_3(t),
    M_p(t) = sum_v rho(p^v) p^{-|v|/2} e^{i v t ln p}: every line r = 2^a 3^b with the law's amplitude (H restricted to
    p = 2, 3: frozen J23)."""
    rng = np.random.default_rng(seed); out = []; clip = []
    for k in range(len(F['ids'])):
        N, w, c = int(F['N'][k]), int(F['w'][k]), float(F['c'][k])
        if kind == 'gram':
            tg = np.linspace(0.5, top, 400001); th = theta_E(N, tg) / math.pi
            phi = 0.5 if w == 1 else 0.0
            kk = np.arange(math.ceil(th[0] - phi) - 1, math.floor(th[-1] - phi) + 2); x0 = kk + phi
            tk = np.interp(x0, th, tg); L = np.log(N * np.maximum(tk, 1.0) ** 2 / (4 * math.pi ** 2))
            j = np.sqrt(np.maximum(np.log(np.maximum(L / c, 1.0001)), 0) / (2 * math.pi ** 2))
            x = x0 + j * rng.normal(size=len(x0)); g = np.interp(x, th, tg)
            out.append(np.sort(g[(g > 0.5) & (g <= top)])); continue
        h = 0.002; t = np.arange(0.5, top, h) + h / 2
        Lt = math.log(N) - 2 * math.log(2 * math.pi) + 2 * np.log(np.abs(1 + 1j * t))       # 2 theta'(t), t >> 1
        dens = Lt / (2 * math.pi)
        if kind == 'lines':
            M = np.ones(len(t), complex)
            for p in (2, 3):
                v = F['M%d_v' % p]; co = F['M%d' % p][k]
                M *= _msum(t, v, co, p)
            a = -w * c * M
            dens = dens + (1 / math.pi) * np.real(a * np.exp(-2j * theta_E(N, t)))
            clip.append(float(np.mean(dens < 0)))
            dens = np.maximum(dens, 0)
        n = rng.poisson(dens * h)
        g = np.repeat(t - h / 2, n) + h * rng.uniform(size=int(n.sum()))
        out.append(np.sort(g))
    if clip:
        print('  synthetic lines: fraction of t with clipped intensity: mean %.4f, max %.4f' % (np.mean(clip), np.max(clip)))
    return out


def _msum(t, v, co, p):
    out = np.zeros(len(t), complex); lp = math.log(p)
    for vv, cc in zip(v, co):
        if abs(cc) > 1e-12: out += cc * np.exp(1j * vv * lp * t)
    return out


# ------------------------------------------------------------------ one-point sums and their lines
def onepoint(F, zeros):
    """per curve: G at xi = ln(a/b) for the registered lines (sharp, Hann) and at the control frequencies."""
    T1, T2 = BAND; T = T2 - T1; L = [tuple(x) for x in F['lines']]
    Gl = np.zeros((len(zeros), len(L), 2), complex); Gc = np.zeros((len(zeros), len(CONTROL_XI), 2), complex)
    for k, g in enumerate(zeros):
        b = g[(g > T1) & (g <= T2)]; ph = np.exp(2j * theta_E(int(F['N'][k]), b)); tap = np.sin(math.pi * (b - T1) / T) ** 2
        for m, (a, bb) in enumerate(L):
            e = ph * np.exp(-1j * math.log(a / bb) * b); Gl[k, m] = e.sum(), (e * tap).sum()
        for m, xi in enumerate(CONTROL_XI):
            e = ph * np.exp(-1j * xi * b); Gc[k, m] = e.sum(), (e * tap).sum()
    return Gl, Gc


def window_lines(F, zeros):
    """one-point sums over the X window W (sharp), for the measured line products."""
    W1, W2 = WIN; L = [tuple(x) for x in F['lines']]; out = np.zeros((len(zeros), len(L)), complex)
    for k, g in enumerate(zeros):
        b = g[(g > W1) & (g <= W2)]; ph = np.exp(2j * theta_E(int(F['N'][k]), b))
        for m, (a, bb) in enumerate(L):
            out[k, m] = np.sum(ph * np.exp(-1j * math.log(a / bb) * b))
    return out


# ------------------------------------------------------------------ X matrix
def xmatrix(F, zeros, sig):
    W1, W2 = WIN; T = W2 - W1
    s = sig / math.sqrt(2); tg = np.arange(W1, W2 + DT / 2, DT); half = int(6 * s / DT) + 1
    G = np.zeros((len(zeros), len(tg)), complex)
    for r, g in enumerate(zeros):
        b = g[(g > W1 - 6 * s) & (g < W2 + 6 * s)]
        ph = np.exp(2j * theta_E(int(F['N'][r]), b)); c = np.round((b - W1) / DT).astype(int)
        for off in range(-half, half + 1):
            j = c + off; ok = (j >= 0) & (j < len(tg))
            np.add.at(G[r], j[ok], ph[ok] * np.exp(-(tg[j[ok]] - b[ok]) ** 2 / (2 * s * s)))
    X = (G @ G.conj().T) * (DT / (math.sqrt(math.pi) * s * T))
    return X


# ------------------------------------------------------------------ regression with jackknife over curves
class Fit:
    """least squares y = A c with the delete-one-curve jackknife (rows carry two curve ids)."""

    def __init__(self, ui, uj, n):
        self.ui, self.uj, self.n = ui, uj, n

    def run(self, A, y, mask=None):
        A = np.asarray(A, float); y = np.asarray(y, float)
        ui, uj = self.ui, self.uj
        if mask is not None: A, y, ui, uj = A[mask], y[mask], ui[mask], uj[mask]
        k = A.shape[1]
        M = A[:, :, None] * A[:, None, :]; v = A * y[:, None]
        MM = M.sum(0); vv = v.sum(0); c = np.linalg.solve(MM, vv)
        SM = np.zeros((self.n, k, k)); Sv = np.zeros((self.n, k))
        np.add.at(SM, ui, M); np.add.at(Sv, ui, v); o = uj != ui
        np.add.at(SM, uj[o], M[o]); np.add.at(Sv, uj[o], v[o])
        present = np.zeros(self.n, bool); present[ui] = True; present[uj] = True
        reps = np.array([np.linalg.solve(MM - SM[f], vv - Sv[f]) for f in np.nonzero(present)[0]])
        nj = len(reps); se = np.sqrt((nj - 1) / nj * ((reps - reps.mean(0)) ** 2).sum(0))
        r = y - A @ c
        return dict(coef=c, se=se, n_rows=int(len(y)), n_units=int(nj), resid_sd=float(np.sqrt(np.mean(r ** 2))))


def readings_pairs(F, X, Gw):
    res = {}; I, J = F['I'], F['J']; n = len(F['ids'])
    fit = Fit(I, J, n); W1, W2 = WIN; T = W2 - W1
    w = F['w'].astype(float); N = F['N'].astype(float); L = [tuple(x) for x in F['lines']]
    ww = w[I] * w[J]; dl = np.abs(np.log(N[J] / N[I]))
    for sig in SIGMAS:
        x = X[sig][I, J]; one = np.ones(len(I)); tg = 's%.1f' % sig
        p = F['pred_%.1f' % sig]; p1 = F['pred1_%.1f' % sig]; pn = F['prednoroot_%.1f' % sig]; pr = F['predrev_%.1f' % sig]
        res[tg + ' X1 slope of Re X on the law'] = fit.run(np.column_stack([one, p]), x.real)
        res[tg + ' X1b slope of Re X on the law, no intercept'] = fit.run(p[:, None], x.real)
        res[tg + ' X2 split [H=1 part, H-1 part]'] = fit.run(np.column_stack([one, p1, p - p1]), x.real)
        res[tg + ' X3 [H=1, H-1 (law), H-1 reversed pairing]'] = fit.run(np.column_stack([one, p1, p - p1, pr - p1]), x.real)
        res[tg + ' X4 [law, no root numbers]'] = fit.run(np.column_stack([one, p, pn]), x.real)
        # measured r = 1 products removed: the blind part
        P0 = sig * math.sqrt(2 * math.pi) * np.real(Gw[I, 0] * np.conj(Gw[J, 0])) / T ** 2
        res[tg + ' X5 Re X - measured (r=1) product, on H-1 part'] = fit.run(np.column_stack([one, p - p1]), x.real - P0)
        # continuum residual: measured products of the registered lines replace their law products
        meas = np.zeros(len(I)); lawp = np.zeros(len(I))
        for m, (a, b) in enumerate(L):
            e = math.exp(-sig * sig * math.log(a / b) ** 2 / 2)
            meas += sig * math.sqrt(2 * math.pi) * e * np.real(Gw[I, m] * np.conj(Gw[J, m])) / T ** 2
            lawp += sig * math.sqrt(2 * math.pi) * e * F['lineamp'][I, m] * F['lineamp'][J, m] / (4 * math.pi ** 2)
        cont = x.real - meas - (p - lawp)
        res[tg + ' X6 continuum: Re X - measured lines - law of the other lines: mean'] = fit.run(one[:, None], cont)
        res[tg + ' X6b continuum: slope on the law'] = fit.run(np.column_stack([one, p]), cont)
        res[tg + ' X7 Im X: mean'] = fit.run(one[:, None], x.imag)
        res[tg + ' X7b Im X: slope on the law'] = fit.run(np.column_stack([one, p]), x.imag)
        res[tg + ' X8 mean Re X (H0: 0)'] = fit.run(one[:, None], x.real)
        res[tg + ' X8b mean of the law (no error)'] = dict(coef=np.array([p.mean()]), se=np.zeros(1), n_rows=len(p), n_units=n, resid_sd=0.0)
        for nm, mm in (('w1 w2 = +1', ww > 0), ('w1 w2 = -1', ww < 0), ('|ln N2/N1| < 0.3', dl < 0.3),
                       ('0.3 <= |ln N2/N1| < 0.8', (dl >= 0.3) & (dl < 0.8)), ('|ln N2/N1| >= 0.8', dl >= 0.8)):
            if mm.sum() > 30:
                res[tg + ' X9 slope, ' + nm] = fit.run(np.column_stack([one, p]), x.real, mm)
        neg = p < -0.02
        if neg.sum() > 30:
            res[tg + ' X10 pairs with law < -0.02: slope, no intercept'] = fit.run(p[:, None], x.real, neg)
    return res


def readings_onepoint(F, Gl, Gc):
    res = {}; n = len(F['ids']); u = np.arange(n); fit = Fit(u, u, n)
    T = BAND[1] - BAND[0]; L = [tuple(x) for x in F['lines']]
    law = F['lineamp'] * T / (2 * math.pi); lawn = F['lineamp_noroot'] * T / (2 * math.pi)
    res['O1 r = 1: slope of Re G (sharp) on -w c T/2pi [not blind]'] = fit.run(law[:, :1], Gl[:, 0, 0].real)
    res['O1b r = 1: slope of Re G (Hann) on the law/2 [not blind]'] = fit.run(law[:, :1] / 2, Gl[:, 0, 1].real)
    m = np.arange(1, len(L))
    y = Gl[:, m, 0].real.ravel(); A = law[:, m].ravel()[:, None]; uu = np.repeat(u, len(m))
    f2 = Fit(uu, uu, n)
    res['O2 lines r != 1 pooled (12 lines): slope of Re G (sharp) on the law'] = f2.run(A, y)
    res['O2b lines r != 1 pooled: slope of Re G (Hann) on the law/2'] = f2.run(A / 2, Gl[:, m, 1].real.ravel())
    res['O2c lines r != 1: [law, law without root number]'] = f2.run(np.column_stack([A[:, 0], lawn[:, m].ravel()]), y)
    for k in m:
        res['O3 line r = %d/%d: slope (sharp)' % L[k]] = fit.run(law[:, k:k + 1], Gl[:, k, 0].real)
    res['O4 Im G, all lines incl. r = 1: mean'] = Fit(np.repeat(u, len(L)), np.repeat(u, len(L)), n).run(np.ones((n * len(L), 1)), Gl[:, :, 0].imag.ravel())
    cc = Gc[:, :, 0].real.ravel(); uc = np.repeat(u, len(CONTROL_XI))
    res['O5 control frequencies: mean Re G'] = Fit(uc, uc, n).run(np.ones((len(cc), 1)), cc)
    res['O5b control frequencies: rms Re G (no error)'] = dict(coef=np.array([np.sqrt(np.mean(cc ** 2))]), se=np.zeros(1), n_rows=len(cc), n_units=n, resid_sd=0.0)
    return res


def fmt(name, r):
    c = np.atleast_1d(r['coef'])
    return '  %-66s %s   [rows %d, units %d]' % (name, ', '.join('%.4f +- %.4f' % (a, b) for a, b in zip(c, r['se'])), r['n_rows'], r['n_units'])


def jsonable(r):
    return {k: (np.atleast_1d(v).tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('step'); ap.add_argument('--design', default='fresh')
    ap.add_argument('--zeros', default=None); ap.add_argument('--frozen', default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument('--out', default=None); ap.add_argument('--synthetic', default=None); ap.add_argument('--seed', type=int, default=1)
    A = ap.parse_args()
    out = A.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_%s_%s' % (A.design, A.synthetic or 'data'))
    os.makedirs(out, exist_ok=True)
    F = load_frozen(A.frozen, A.design); t0 = time.time()
    steps = ['onepoint', 'x', 'readings'] if A.step == 'all' else [A.step]
    zeros = None
    def Z():
        nonlocal zeros
        if zeros is None:
            zeros = synth_zeros(F, A.synthetic, A.seed) if A.synthetic else load_zeros(A.design, F, A.zeros or DEFAULT_ZEROS[A.design])
        return zeros
    fo = os.path.join(out, 'onepoint.npz')
    if 'onepoint' in steps and not os.path.exists(fo):
        Gl, Gc = onepoint(F, Z()); Gw = window_lines(F, Z()); np.savez(fo, Gl=Gl, Gc=Gc, Gw=Gw)
        print('one-point sums done (%.0fs)' % (time.time() - t0), flush=True)
    for sig in SIGMAS:
        fx = os.path.join(out, 'X_%.1f.npy' % sig)
        if 'x' in steps and not os.path.exists(fx):
            np.save(fx, xmatrix(F, Z(), sig)); print('X sigma %.1f done (%.0fs)' % (sig, time.time() - t0), flush=True)
    if 'readings' in steps:
        O = np.load(fo); X = {s: np.load(os.path.join(out, 'X_%.1f.npy' % s)) for s in SIGMAS}
        ro = readings_onepoint(F, O['Gl'], O['Gc']); rp = readings_pairs(F, X, O['Gw'])
        print('readings (%s, %s, %d curves, %d pairs):' % (A.design, A.synthetic or 'data', len(F['ids']), len(F['I'])))
        for k, v in list(ro.items()) + list(rp.items()): print(fmt(k, v))
        json.dump({k: jsonable(v) for k, v in list(ro.items()) + list(rp.items())}, open(os.path.join(out, 'readings.json'), 'w'), indent=1)
    print('done (%.0fs)' % (time.time() - t0))


if __name__ == '__main__':
    main()
