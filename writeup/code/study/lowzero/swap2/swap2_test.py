#!/usr/bin/env python3
"""swap2: registered re-test of the swapped ("Gram-lattice") term between the zeros of two Dirichlet L-functions.

Statistic (per ordered pair of one-sided L-functions i, j; phases cancel the recipe's phase exactly):
   X_ij(sigma) = (1/T) sum_{g_i, g_j} e^{2i theta_i(g_i) - 2i theta_j(g_j)} K(g_j - g_i) w((g_i + g_j)/2),
   computed as (1/(sqrt(pi) s T)) int_W G_i(t) conj G_j(t) dt with G_i(t) = sum_g e^{2i theta_i(g)} exp(-(t - g)^2 / 2 s^2),
   s = sigma/sqrt 2, K(u) = exp(-u^2/2 sigma^2), W = [W1, W2] the design's window (6 units inside the stored band at
   every edge without stored zeros), T = W2 - W1, t-grid step 0.02.
   theta_chi(t) = (t/2) ln(q/pi) + Im lnGamma((1/2 + a + i t)/2)  (a = parity; exact loggamma).
   One-sided functions: real chi_D (zeros gamma > 0); complex chi (zeros gamma > 0) and conj chi (zeros -gamma, gamma < 0).
Law (frozen in pred_*.npz):  E X_ij = (1/4 pi^2) k_i k_j eps_i conj(eps_j) J_ij,  J_ij = int H_ij(u) K(u) du.
One-point (frozen in onepoint.npz):  G_chi = sum_{T1 < g <= T2} e^{2i theta_chi(g)} ~ -eps k (T2 - T1)/2pi; per function
   the resonance prime sum R_w = -(i^a/sqrt q) sum_n chi(n) Lambda(n) e(-n/q) w(2 pi n/q) (sharp and Hann windows).

usage:  swap2_test.py STEP [--zeros DIR] [--frozen DIR] [--out DIR] [--designs OO,HI,C] [--synthetic poisson|gram] [--seed S]
  STEP = complete | onepoint | x | readings | all.  Every part is cached in --out (npy/json) and skipped when present,
  so an interrupted run resumes; each part takes less than 10 minutes on 2 cores.
  Designs: OO (primary; 546 real, 100-600), HI (R120 at 600-1000), C (88 complex characters, both half-lines, and 60 real,
  100-500), FN (VALIDATION ONLY: the used set of round 138, where the retrospective 1.004 was measured).
Needs numpy, scipy; cypari2 only for the optional root-number check (step complete)."""
import sys, os, math, json, time, hashlib, argparse
import numpy as np
from scipy.special import loggamma

FROZEN_SHA = {   # sha256 of the frozen prediction files (filled at freezing; the test refuses to run on other files)
    'pred_OO.npz': '884866caff9ef4af29a264889e93372ef3c00f0d911717006168ed47fba7c4a8',
    'pred_C.npz': 'f1fa89eda1f1760f0cfabc0bfc6374a30c1ccf853ef1171601ffe1da2acd249a',
    'pred_FN.npz': '75098b6ab49936d48b8ed899718228b5767d5c0c5e22ff52db8b7563cf6bac21',
    'onepoint.npz': '2ee4adf440443b4e78be7bf952852952525cb4e2f93374c1f90456c318bd905d',
}
SIGMAS = (0.5, 1.0)
DT = 0.02
ZFILES = {
    'R120': ['r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'],
    'B': ['r134_B.npz'], 'DN': ['r136_Dn.npz'], 'C134': ['r134_C.npz'],
    'CPLX': ['r121_part0.npz', 'r121_part1.npz', 'r121_part2.npz', 'r121_part3.npz'], 'FN': ['r138_Fn.npz'],
}
DESIGN = {   # band = stored heights used for completeness and the one-point sum; win = X window
    'OO': dict(band=(100.0, 600.0), win=(100.0, 594.0)),
    'HI': dict(band=(600.0, 1000.0), win=(600.0, 994.0)),
    'C': dict(band=(100.0, 500.0), win=(106.0, 494.0)),
    'FN': dict(band=(100.0, 600.0), win=(100.0, 594.0)),
}
# completeness thresholds: real sets (count from 0; round 138): |N - theta/pi| <= 1.5, largest unfolded gap <= 3.0;
# one-sided sets stored from 100 only (complex): centred count <= 1.5, step statistic <= 0.7, gap <= 3.0
TH_COUNT, TH_GAP, TH_STEP = 1.5, 3.0, 0.7


def theta(t, q, a):
    t = np.asarray(t, float)
    return t / 2 * math.log(q / math.pi) + np.imag(loggamma(0.25 + a / 2 + 0.5j * t))


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def phi_ratio(k):
    r = 1.0; n = k; p = 2
    while p * p <= n:
        if n % p == 0:
            r *= p / (p - 1)
            while n % p == 0: n //= p
        p += 1
    if n > 1: r *= n / (n - 1)
    return r


# ------------------------------------------------------------------ functions of each design
def load_frozen(fdir):
    F = {}
    for f in ('pred_OO.npz', 'pred_C.npz', 'pred_FN.npz', 'onepoint.npz'):
        p = os.path.join(fdir, f)
        if FROZEN_SHA.get(f) and sha(p) != FROZEN_SHA[f]:
            raise SystemExit('frozen file %s does not match its registered sha256' % f)
        F[f] = dict(np.load(p, allow_pickle=False))
    return F


def design_functions(des, F):
    """list of dicts: label, kind ('R' or 'C'), D or (N, m, conj), q, a, eps, kappa, unit."""
    out = []
    if des in ('OO', 'HI', 'FN'):
        Ds = F['pred_FN.npz' if des == 'FN' else 'pred_OO.npz']['Ds']
        if des == 'HI': Ds = Ds[:181]
        for D in Ds:
            D = int(D); out.append(dict(label='D%d' % D, kind='R', D=D, q=abs(D), a=0 if D > 0 else 1, eps=1 + 0j, kappa=phi_ratio(abs(D)), unit='D%d' % D))
    else:
        C = F['pred_C.npz']
        for lab, q, a, e, k in zip(C['labels'], C['q'], C['a'], C['eps'], C['kappa']):
            lab = str(lab)
            if lab.startswith('D'):
                out.append(dict(label=lab, kind='R', D=int(lab[1:]), q=int(q), a=int(a), eps=complex(e), kappa=float(k), unit=lab))
            else:
                p_ = lab.split('_'); N, m = int(p_[0][1:]), int(p_[1][1:]); cj = lab.endswith('_conj')
                out.append(dict(label=lab, kind='C', N=N, m=m, conj=cj, q=int(q), a=int(a), eps=complex(e), kappa=float(k), unit='N%d_m%d' % (N, m)))
    return out


def load_zeros(des, funcs, zdir, synthetic=None, seed=0):
    """one-sided zero arrays (sorted, > 0) for every function of the design."""
    if synthetic:
        return synth_zeros(des, funcs, synthetic, seed)
    Z = {}
    groups = {'OO': ['R120', 'B', 'DN'], 'HI': ['C134'], 'FN': ['FN'], 'C': ['CPLX', 'R120']}[des]
    for g in groups:
        for f in ZFILES[g]:
            with np.load(os.path.join(zdir, f)) as z:
                for k in z.files: Z[k] = np.array(z[k], float)
    out = []
    for fn in funcs:
        if fn['kind'] == 'R':
            out.append(np.sort(Z['D%d' % fn['D']]))
        else:
            z = Z['N%d_m%d' % (fn['N'], fn['m'])]
            out.append(np.sort(-z[z < 0]) if fn['conj'] else np.sort(z[z > 0]))
    return out


def synth_zeros(des, funcs, kind, seed):
    """synthetic one-sided sets with each function's exact mean density: 'poisson' (independent points, uniform in
    theta/pi) or 'gram' (one point per Gram interval at theta - arg(eps)/2 = pi (k + 1/2) + pi j xi, xi ~ N(0, 1), with
    e^{-2 pi^2 j^2} = kappa/L(t): reproduces the one-point law exactly and the H = 1 version of the swapped term)."""
    rng = np.random.default_rng(seed)
    lo, hi = DESIGN[des]['band']
    out = []
    for fn in funcs:
        q, a, e, kap = fn['q'], fn['a'], fn['eps'], fn['kappa']
        t0 = 100.0 if (des == 'C' and fn['kind'] == 'C') else 1.0
        tg = np.linspace(t0, hi + 2, 200000); th = theta(tg, q, a) / math.pi
        if kind == 'poisson':
            n = rng.poisson(th[-1] - th[0]); x = np.sort(rng.uniform(th[0], th[-1], n))
        else:
            ph = math.atan2(e.imag, e.real) / (2 * math.pi)
            k = np.arange(math.ceil(th[0] - 0.5 - ph) - 1, math.floor(th[-1] - 0.5 - ph) + 2)
            x0 = k + 0.5 + ph
            tk = np.interp(x0, th, tg); L = np.log(q * np.maximum(tk, 1.0) / (2 * math.pi))
            j = np.sqrt(np.maximum(np.log(np.maximum(L / kap, 1.0001)), 0) / (2 * math.pi ** 2))
            x = x0 + j * rng.normal(size=len(x0))
        g = np.interp(x, th, tg)
        out.append(np.sort(g[(g > t0) & (g <= hi + 2)]))
    return out


# ------------------------------------------------------------------ completeness
def completeness(des, funcs, zeros):
    rows = []
    for fn, g in zip(funcs, zeros):
        q, a = fn['q'], fn['a']
        if des == 'C' and fn['kind'] == 'R':
            lo, hi = 100.0, 600.0
        else:
            lo, hi = DESIGN[des]['band']
        t = np.linspace(lo, hi, 200); b = g[(g >= lo) & (g <= hi)]
        u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e)); gap = float(np.diff(u).max())
        if fn['kind'] == 'R':
            s = np.searchsorted(g, t) - theta(t, q, a) / math.pi
            cnt = float(np.abs(s).max()); step = float('nan')
            flag = cnt > TH_COUNT or gap > TH_GAP
        else:
            s = np.searchsorted(b, t, 'right') - (theta(t, q, a) - theta(lo, q, a)) / math.pi
            cnt = float(np.abs(s - s.mean()).max()); cs = np.cumsum(s)
            step = float(max(abs(cs[k - 1] / k - (cs[-1] - cs[k - 1]) / (200 - k)) for k in range(20, 181)))
            flag = cnt > TH_COUNT or step > TH_STEP or gap > TH_GAP
        rows.append(dict(label=fn['label'], count=cnt, step=step, gap=gap, n_band=int(len(b)), flagged=bool(flag)))
    return rows


# ------------------------------------------------------------------ one-point sums
def onepoint(des, funcs, zeros):
    T1, T2 = DESIGN[des]['band']
    out = []
    for fn, g in zip(funcs, zeros):
        b = g[(g > T1) & (g <= T2)]; ph = np.exp(2j * theta(b, fn['q'], fn['a']))
        x = (b - T1) / (T2 - T1)
        out.append((complex(ph.sum()), complex((ph * np.sin(math.pi * x) ** 2).sum())))
    return np.array(out)


# ------------------------------------------------------------------ X matrix
def xmatrix(des, funcs, zeros, sig, rows=None):
    W1, W2 = DESIGN[des]['win']; T = W2 - W1
    s = sig / math.sqrt(2); tg = np.arange(W1, W2 + DT / 2, DT); half = int(6 * s / DT) + 1
    idx = range(len(funcs)) if rows is None else rows
    G = np.zeros((len(idx), len(tg)), complex)
    for r, i in enumerate(idx):
        fn, g = funcs[i], zeros[i]
        b = g[(g > W1 - 6 * s) & (g < W2 + 6 * s)]
        ph = np.exp(2j * theta(b, fn['q'], fn['a'])); c = np.round((b - W1) / DT).astype(int)
        for off in range(-half, half + 1):
            j = c + off; ok = (j >= 0) & (j < len(tg))
            np.add.at(G[r], j[ok], ph[ok] * np.exp(-(tg[j[ok]] - b[ok]) ** 2 / (2 * s * s)))
    return G, DT / (math.sqrt(math.pi) * s * T)


# ------------------------------------------------------------------ regression with jackknife over units
class Fit:
    """least squares y = A c (complex allowed) with the delete-one-unit jackknife; rows carry two unit ids."""

    def __init__(self, ui, uj, nunits):
        self.ui, self.uj, self.n = ui, uj, nunits

    def run(self, A, y, mask=None):
        A = np.asarray(A); y = np.asarray(y)
        if mask is not None: A, y, ui, uj = A[mask], y[mask], self.ui[mask], self.uj[mask]
        else: ui, uj = self.ui, self.uj
        cplx = np.iscomplexobj(A) or np.iscomplexobj(y)
        dt = complex if cplx else float
        A = A.astype(dt); y = y.astype(dt); k = A.shape[1]
        M = A.conj()[:, :, None] * A[:, None, :]; v = A.conj() * y[:, None]
        MM = M.sum(0); vv = v.sum(0); c = np.linalg.solve(MM, vv)
        SM = np.zeros((self.n, k, k), dt); Sv = np.zeros((self.n, k), dt)
        np.add.at(SM, ui, M); np.add.at(Sv, ui, v)
        o = uj != ui
        np.add.at(SM, uj[o], M[o]); np.add.at(Sv, uj[o], v[o])
        present = np.zeros(self.n, bool); present[ui] = True; present[uj] = True
        reps = np.array([np.linalg.solve(MM - SM[f], vv - Sv[f]) for f in np.nonzero(present)[0]])
        nj = len(reps)
        se_re = np.sqrt((nj - 1) / nj * ((reps.real - reps.real.mean(0)) ** 2).sum(0))
        se_im = np.sqrt((nj - 1) / nj * ((reps.imag - reps.imag.mean(0)) ** 2).sum(0)) if cplx else np.zeros(k)
        r = y - A @ c
        return dict(coef=c, se_re=se_re, se_im=se_im, n_rows=int(len(y)), n_units=int(nj),
                    resid_sd=float(np.sqrt(np.mean(np.abs(r) ** 2))))


def readings_pairs(des, F, X, excl, funcs):
    """all registered pair readings for one design and one sigma-indexed X dict; returns {name: result}."""
    res = {}
    if des in ('OO', 'HI', 'FN'):
        P = F['pred_FN.npz' if des == 'FN' else 'pred_OO.npz']
        I, J = P['I'], P['J']
        if des == 'HI':
            m = (I < 181) & (J < 181); I, J = I[m], J[m]
        else:
            m = np.ones(len(I), bool)
        keep = ~(excl[I] | excl[J])
        units = {f['unit']: k for k, f in enumerate(funcs)}; uid = np.array([units[f['unit']] for f in funcs])
        fit = Fit(uid[I][keep], uid[J][keep], len(units))
        q = np.array([f['q'] for f in funcs]); a = np.array([f['a'] for f in funcs])
        dl = np.abs(np.log(q[J] / q[I]))[keep]; par = (a[I] + a[J])[keep]
        for sig in SIGMAS:
            x = X[sig][I, J][keep]; p = P['pred_%.1f' % sig][m][keep]; p1 = P['pred1_%.1f' % sig][m][keep]; one = np.ones(len(x))
            tag = 's%.1f' % sig
            res[tag + ' R1 slope of Re X on the law'] = fit.run(np.column_stack([one, p]), x.real)
            res[tag + ' R2 split [H=1 part, H-1 part]'] = fit.run(np.column_stack([one, p1, p - p1]), x.real)
            for lo, hi in ((0, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 9.0)):
                mm = (dl >= lo) & (dl < hi)
                if mm.sum() > 50:
                    res[tag + ' R3 slope, %.1f <= |ln(q2/q1)| < %.1f' % (lo, hi)] = fit.run(np.column_stack([one, p]), x.real, mm)
            res[tag + ' R4 slope of Im X on the law'] = fit.run(np.column_stack([one, p]), x.imag)
            res[tag + ' R4b mean Im X'] = fit.run(one[:, None], x.imag)
            for nm, mm in (('even-even', par == 0), ('odd-odd', par == 2), ('even-odd', par == 1)):
                res[tag + ' R5 slope, ' + nm] = fit.run(np.column_stack([one, p]), x.real, mm)
            neg = p < -0.01
            res[tag + ' R6 pairs with law < -0.01: slope of Re X on the law, no intercept'] = fit.run(p[:, None], x.real, neg)
            res[tag + ' R7 slope of Re X on the law, no intercept'] = fit.run(p[:, None], x.real)
    else:
        P = F['pred_C.npz']
        units = {}; [units.setdefault(f['unit'], len(units)) for f in funcs]; uid = np.array([units[f['unit']] for f in funcs])
        for sig in SIGMAS:
            tag = 's%.1f' % sig
            blocks = []
            for nm in ('cc', 'cr'):
                I, J = P['I_%s' % nm], P['J_%s' % nm]
                blocks.append((nm, I, J, P['pred_%s_%.1f' % (nm, sig)], P['pred1_%s_%.1f' % (nm, sig)], P['prednoroot_%s_%.1f' % (nm, sig)], P['J_%s_%.1f' % (nm, sig)]))
            I = np.concatenate([b[1] for b in blocks]); J = np.concatenate([b[2] for b in blocks])
            p = np.concatenate([b[3] for b in blocks]); p1 = np.concatenate([b[4] for b in blocks]); pn = np.concatenate([b[5] for b in blocks])
            JJ = np.concatenate([b[6] for b in blocks]); iscc = np.concatenate([np.full(len(b[1]), b[0] == 'cc') for b in blocks])
            keep = ~(excl[I] | excl[J]); I, J, p, p1, pn, JJ, iscc = I[keep], J[keep], p[keep], p1[keep], pn[keep], JJ[keep], iscc[keep]
            x = X[sig][I, J]; one = np.ones(len(x), complex)
            kk = np.array([f['kappa'] for f in funcs]); ee = np.array([f['eps'] for f in funcs])
            base = kk[I] * kk[J] * ee[I] * np.conj(ee[J]) / (4 * math.pi ** 2)
            fit = Fit(uid[I], uid[J], len(units))
            mirror = np.array([funcs[i]['unit'] == funcs[j]['unit'] for i, j in zip(I, J)])
            res[tag + ' C1 complex slope of X on the law (CC + CR)'] = fit.run(np.column_stack([one, p]), x)
            res[tag + ' C2 two predictors [law, no root numbers]'] = fit.run(np.column_stack([one, p, pn]), x)
            res[tag + ' C3 [Re J part, i Im J part] (direction: law 1, conj J -1)'] = fit.run(np.column_stack([one, base * JJ.real, base * 1j * JJ.imag]), x)
            res[tag + ' C4 (chi, conj chi) pairs: complex slope'] = fit.run(np.column_stack([one, p]), x, mirror)
            res[tag + ' C5 complex x real pairs: complex slope'] = fit.run(np.column_stack([one, p]), x, ~iscc)
            res[tag + ' C6 complex x complex, chi2 != conj chi1: complex slope'] = fit.run(np.column_stack([one, p]), x, iscc & ~mirror)
            res[tag + ' C7 split [H=1 part, H-1 part] (CC + CR)'] = fit.run(np.column_stack([one, p1, p - p1]), x)
    return res


def readings_onepoint(des, F, G, excl, funcs):
    O = F['onepoint.npz']; key = 'FN' if des == 'FN' else ('HI' if des == 'HI' else ('C' if des == 'C' else 'OO'))
    labs = [str(s) for s in O[key + '_labels']]; idx = {l: k for k, l in enumerate(labs)}
    sel = [k for k, f in enumerate(funcs) if f['label'] in idx and not excl[k]]
    k_ = np.array([idx[funcs[k]['label']] for k in sel])
    Gs, Gh = G[sel, 0], G[sel, 1]
    Ls, Lh, Rs, Rh = O[key + '_law_sharp'][k_], O[key + '_law_hann'][k_], O[key + '_res_sharp'][k_], O[key + '_res_hann'][k_]
    units = {}; [units.setdefault(funcs[k]['unit'], len(units)) for k in sel]; uid = np.array([units[funcs[k]['unit']] for k in sel])
    fit = Fit(uid, uid, len(units)); res = {}
    a = np.array([funcs[k]['a'] for k in sel]); kap = np.array([funcs[k]['kappa'] for k in sel])
    T1, T2 = DESIGN[des]['band']
    if des == 'C':
        isc = np.array([funcs[k]['kind'] == 'C' for k in sel])
        Gs, Gh, Ls, Lh, Rs, Rh, a, kap, uid = [v[isc] for v in (Gs, Gh, Ls, Lh, Rs, Rh, a, kap, uid)]
        fit = Fit(uid, uid, len(units))
        res['P4 complex slope of G (sharp) on -eps k T/2pi'] = fit.run(Ls[:, None], Gs)
        res['P4b complex slope of G (sharp) on -k T/2pi (no root number)'] = fit.run(np.abs(Ls)[:, None] * -1, Gs)
        res['P4c prime side: complex slope of R (sharp) on the law [frozen value, no error]'] = fit.run(Ls[:, None], Rs)
        res['P5 data check: complex slope of G (Hann) on R (Hann)'] = fit.run(Rh[:, None], Gh)
    else:
        res['P1 slope of Re G (sharp) on -k T/2pi'] = fit.run(Ls[:, None].real, Gs.real)
        res['P1b prime side: slope of Re R (sharp) on the law [frozen value]'] = fit.run(Ls[:, None].real, Rs.real)
        res['P1c slope of Re G (Hann) on the law'] = fit.run(Lh[:, None].real, Gh.real)
        res['P1d prime side: slope of Re R (Hann) on the law [frozen value]'] = fit.run(Lh[:, None].real, Rh.real)
        res['P2 data check: slope of Re G (Hann) on Re R (Hann)'] = fit.run(Rh[:, None].real, Gh.real)
        res['P3 even: slope of Re G (sharp) on the law'] = fit.run(Ls[:, None].real, Gs.real, a == 0)
        res['P3 odd: slope of Re G (sharp) on the law'] = fit.run(Ls[:, None].real, Gs.real, a == 1)
        y = Gs.real * 2 * math.pi / (T2 - T1)
        res['P3b Re G 2pi/T = c0 + c1 kappa (law: 0, -1)'] = fit.run(np.column_stack([np.ones(len(y)), kap]), y)
    d = np.abs(Gh - Rh)
    res['P6 |G - R| (Hann): sd, max, count > 1'] = dict(coef=np.array([np.sqrt(np.mean(d ** 2)), d.max(), float((d > 1).sum())]),
                                                     se_re=np.zeros(3), se_im=np.zeros(3), n_rows=int(len(d)), n_units=len(units), resid_sd=0.0)
    return res


def fmt(name, r):
    c = np.atleast_1d(r['coef'])
    parts = []
    for v, e1, e2 in zip(c, r['se_re'], r['se_im']):
        if abs(complex(v).imag) > 0 or e2 > 0:
            parts.append('%.4f%+.4fi (+- %.4f, %.4f)' % (complex(v).real, complex(v).imag, e1, e2))
        else:
            parts.append('%.4f +- %.4f' % (complex(v).real, e1))
    return '  %-62s %s   [rows %d, units %d]' % (name, ', '.join(parts), r['n_rows'], r['n_units'])


def jsonable(r):
    return {k: ([[float(complex(x).real), float(complex(x).imag)] for x in np.atleast_1d(v)] if k == 'coef' else (v.tolist() if isinstance(v, np.ndarray) else v)) for k, v in r.items()}


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('step'); ap.add_argument('--zeros', default=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')))
    ap.add_argument('--frozen', default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument('--out', default=None); ap.add_argument('--designs', default='OO,HI,C')
    ap.add_argument('--synthetic', default=None); ap.add_argument('--seed', type=int, default=1)
    A = ap.parse_args()
    out = A.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_' + (A.synthetic or 'data'))
    os.makedirs(out, exist_ok=True)
    F = load_frozen(A.frozen)
    t0 = time.time()
    steps = ['complete', 'onepoint', 'x', 'readings'] if A.step == 'all' else [A.step]
    report = {}
    for des in A.designs.split(','):
        funcs = design_functions(des, F)
        zeros = None
        def Z():
            nonlocal zeros
            if zeros is None: zeros = load_zeros(des, funcs, A.zeros, A.synthetic, 1000 * A.seed + {'OO': 1, 'HI': 2, 'C': 3, 'FN': 4}[des])
            return zeros
        fc = os.path.join(out, 'complete_%s.json' % des)
        if 'complete' in steps and not os.path.exists(fc):
            rows = completeness(des, funcs, Z())
            json.dump(rows, open(fc, 'w'), indent=0)
            if des == 'C':
                try:
                    import cypari2
                    pari = cypari2.Pari(); pari.allocatemem(200 * 10 ** 6)
                    d = max(abs(complex(pari('lfunrootres(lfuncreate(Mod(%d,%d)))[3]' % (f['m'], f['N']))) - (f['eps'].conjugate() if f['conj'] else f['eps']))
                            for f in funcs if f['kind'] == 'C')
                    print('  root numbers: max |eps(PARI) - frozen eps| = %.1e' % d)
                except ImportError:
                    print('  (cypari2 not available: root numbers not rechecked)')
        if not os.path.exists(fc): raise SystemExit('run step complete first')
        rows = json.load(open(fc))
        excl = np.array([r['flagged'] for r in rows])
        if A.synthetic: excl[:] = False          # synthetic sets are not zeros: flags reported, nothing excluded
        flagged = [(r['label'], round(r['count'], 2), round(r['step'], 2) if r['step'] == r['step'] else None, round(r['gap'], 2)) for r in rows if r['flagged']]
        print('%s: %d one-sided sets, %d flagged%s: %s' % (des, len(rows), len(flagged), ' (synthetic: not excluded)' if A.synthetic else '', flagged[:12]), flush=True)
        fo = os.path.join(out, 'onepoint_%s.npy' % des)
        if 'onepoint' in steps and not os.path.exists(fo):
            np.save(fo, onepoint(des, funcs, Z()))
        for sig in SIGMAS:
            fx = os.path.join(out, 'X_%s_%.1f.npy' % (des, sig))
            if 'x' in steps and not os.path.exists(fx):
                G, scale = xmatrix(des, funcs, Z(), sig)
                X = np.zeros((len(funcs), len(funcs)), complex); Gh = G.conj().T
                for r0 in range(0, len(funcs), 64):
                    X[r0:r0 + 64] = (G[r0:r0 + 64] @ Gh) * scale
                del G, Gh
                np.save(fx, X)
                print('  X %s sigma %.1f done (%.0fs)' % (des, sig, time.time() - t0), flush=True)
        if 'readings' in steps:
            X = {s: np.load(os.path.join(out, 'X_%s_%.1f.npy' % (des, s))) for s in SIGMAS}
            G = np.load(fo)
            rp = readings_pairs(des, F, X, excl, funcs); ro = readings_onepoint(des, F, G, excl, funcs)
            print('%s readings (%s):' % (des, A.synthetic or 'data'))
            for k, v in list(rp.items()) + list(ro.items()): print(fmt(k, v))
            report[des] = {k: jsonable(v) for k, v in list(rp.items()) + list(ro.items())}
            report[des + '_excluded'] = flagged
    if report:
        json.dump(report, open(os.path.join(out, 'readings.json'), 'w'), indent=1)
    print('done (%.0fs)' % (time.time() - t0))


if __name__ == '__main__':
    main()
