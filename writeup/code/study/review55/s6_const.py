"""Own constants c_L = prod_p G_p(1/p), G_p = F_p/(1 - X), F_p(X) = sum_k conj(lam(p^k)) mu(p^k) X^k, from local polynomials:
  T2 newforms: a_p from mfcoefs (good: 1 - lam X + X^2; bad p || N: 1 - lam_p X, lam_p = a_p / p^((k-1)/2));
  T3/T4a Sym^m E (semistable): a_p(E) from ellan; good: prod_j (1 - alpha^(m-2j) X); multiplicative: 1 - (a_p/sqrt p)^m X;
  T4b genus 2: good p: hyperellcharpoly; bad p: PARI lfuneuler (the L-function that produced the zeros).
Partial products at checkpoints; tail estimate nu/(sqrt(P) ln P) (prime-square bias of the adjoint, nu = number of orthogonal
constituents of Ad); log-Cesaro mean over [P/10, P]."""
import json, math, os, sys, time
import numpy as np, cypari2
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'; D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json'))
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)

def F_of(alpha_rows, X):
    """alpha_rows: (n, d) complex Satake parameters (rows may contain zeros for bad p); X: (n,) -> F(X) (n,) complex.
    F(X) = sum_k conj(h_k(alpha)) (-1)^k e_k(alpha) X^k, k = 0..d."""
    n, d = alpha_rows.shape
    # e_k via polynomial prod (1 - a X): coefficients c_k = (-1)^k e_k
    c = np.zeros((n, d + 1), complex); c[:, 0] = 1
    for j in range(d):
        a = alpha_rows[:, j]
        c[:, 1:] = c[:, 1:] - a[:, None] * c[:, :-1]
    # h_k: series of 1/prod(1 - a X): h_0 = 1, h_k = -sum_{i=1..k} c_i h_{k-i}
    h = np.zeros((n, d + 1), complex); h[:, 0] = 1
    for k in range(1, d + 1):
        h[:, k] = -sum(c[:, i] * h[:, k - i] for i in range(1, k + 1))
    return sum(np.conj(h[:, k]) * c[:, k] * X ** k for k in range(d + 1))

def satake_sym(ap, p, m):
    """good p: alpha = e^{i phi}, cos phi = a_p / (2 sqrt p)."""
    x = np.clip(ap / (2 * np.sqrt(p)), -1, 1); phi = np.arccos(x)
    return np.stack([np.exp(1j * (m - 2 * j) * phi) for j in range(m + 1)], axis=1)

def report(key, ps, logG, nu, c_reg, extra=''):
    cum = np.cumsum(logG)
    out = {}
    for Pc in (1e3, 1e4, 1e5, 3e5, 1e6):
        i = np.searchsorted(ps, Pc, side='right') - 1
        if ps[-1] < Pc * 0.99: continue
        cP = math.exp(cum[i]); corr = math.exp(cum[i] + nu / (math.sqrt(Pc) * math.log(Pc)))
        # log-Cesaro mean over [P/10, P] of the corrected log partial products
        sel = (ps >= Pc / 10) & (ps <= Pc)
        lc = cum[sel] + nu / (np.sqrt(ps[sel]) * np.log(ps[sel]))
        wts = 1 / ps[sel]
        ces = math.exp(np.sum(wts * lc) / np.sum(wts))
        out[int(Pc)] = (cP, corr, ces)
    Pm = max(out)
    print('%-28s reg c = %8.5f | ' % (key, c_reg) + '  '.join('P=%.0e: raw %.5f corr %.5f ces %.5f' % (k, *v) for k, v in sorted(out.items()) if k >= 1e4)
          + '  | best/reg = %.4f %s' % (out[Pm][2] / c_reg if c_reg else float('nan'), extra), flush=True)
    return out

def const_sym(key, e, m, Pmax):
    a = e['spec'][1]; E = pari.ellinit(a); N = int(pari.ellglobalred(E)[0])
    an = np.array([int(x) for x in pari.ellan(E, Pmax)], dtype=float)
    ps = np.array([int(x) for x in pari('primes([2,%d])' % Pmax)], dtype=np.int64)
    ap = an[ps - 1]
    good = (N % ps) != 0
    alpha = np.zeros((len(ps), m + 1), complex)
    alpha[good] = satake_sym(ap[good], ps[good].astype(float), m)
    bad = ~good
    alpha[bad, 0] = (ap[bad] / np.sqrt(ps[bad].astype(float))) ** m   # multiplicative reduction only (semistable)
    if np.any(np.abs(ap[bad]) != 1): print('  WARNING additive prime in', key, ps[bad][np.abs(ap[bad]) != 1])
    X = 1 / ps.astype(float)
    G = F_of(alpha, X) / (1 - X)
    nu = {2: 2, 3: 3}[m]           # Ad Sym^2 = Sym^2 + Sym^4 ; Ad Sym^3 = Sym^2 + Sym^4 + Sym^6 (all orthogonal)
    # cross-check bad-prime factor against PARI lfuneuler for the first bad prime
    L = pari.lfunsympow(E, m); pb = int(ps[bad][0]); eul = pari.lfuneuler(L, pb)
    return ps, np.log(G.real), nu, float(np.max(np.abs(G.imag))), (pb, str(eul), alpha[bad][0, 0])

def const_mf(key, e, Pmax):
    k, N, idx = e['spec'][1], e['spec'][2], e['spec'][3]
    mf = pari.mfinit([N, k], 0); F = pari.mfeigenbasis(mf)[idx]
    co = pari.lfunan(pari.lfunmf(mf, F), Pmax)
    ps = np.array([int(x) for x in pari('primes([2,%d])' % Pmax)], dtype=np.int64)
    ap = np.array([float(co[int(p) - 1]) for p in ps])
    lam = ap / ps.astype(float) ** ((k - 1) / 2)
    good = (N % ps) != 0
    alpha = np.zeros((len(ps), 2), complex)
    phi = np.arccos(np.clip(lam[good] / 2, -1, 1)); alpha[good, 0] = np.exp(1j * phi); alpha[good, 1] = np.exp(-1j * phi)
    alpha[~good, 0] = lam[~good]
    X = 1 / ps.astype(float); G = F_of(alpha, X) / (1 - X)
    return ps, np.log(G.real), 1, float(np.max(np.abs(G.imag))), None

def const_g2(key, e, Pmax):
    f = pari(e['spec'][1]); h = pari(e['spec'][2]); L = pari('lfungenus2([%s,%s])' % (e['spec'][1], e['spec'][2]))
    N = int(e['Q'])
    ps = np.array([int(x) for x in pari('primes([2,%d])' % Pmax)], dtype=np.int64)
    alpha = np.zeros((len(ps), 4), complex)
    badlist = []
    for i, p in enumerate(ps):
        p = int(p)
        if N % p == 0 or p == 2:
            # local factor from PARI's L-function (also p = 2 to avoid char-2 subtleties): 1/poly in X = p^-s (arith)
            eu = pari.lfuneuler(L, p); den = pari.denominator(eu) if p else None
            pol = [complex(x) for x in pari.Vec(pari.denominator(eu))] if True else None
            # Vec gives coefficients from highest degree; denominator is poly in x = p^{-s}
            co = np.array(pol[::-1])            # c_0 + c_1 x + ...
            co = co / co[0]
            # analytic normalisation: x = X / sqrt(p)  (weight 2, (k-1)/2 = 1/2)
            coA = co * (p ** -0.5) ** np.arange(len(co))
            roots = np.roots(coA[::-1]) if len(coA) > 1 else np.array([])
            al = 1 / roots if len(roots) else np.array([])
            alpha[i, :len(al)] = al
            badlist.append((p, [round(x.real, 6) for x in co]))
        else:
            cp = pari.hyperellcharpoly([f * pari.Mod(1, p), h * pari.Mod(1, p)])
            v = [float(x) for x in pari.Vec(cp)]          # x^4 - a1 x^3 + a2 x^2 - p a1 x + p^2
            r = np.roots(v) / math.sqrt(p)                  # Frobenius eigenvalues / sqrt p = Satake parameters
            alpha[i] = r
    X = 1 / ps.astype(float); G = F_of(alpha, X) / (1 - X)
    return ps, np.log(G.real), 2, float(np.max(np.abs(G.imag))), badlist

if __name__ == '__main__':
    which = sys.argv[1]; Pmax = int(float(sys.argv[2])); keys = sys.argv[3:]
    res = {}
    for key in keys:
        e = P['L'][key]; t0 = time.time()
        if which == 'mf': ps, lg, nu, im, x = const_mf(key, e, Pmax)
        elif which in ('sym2', 'sym3'): ps, lg, nu, im, x = const_sym(key, e, int(which[-1]), Pmax)
        else: ps, lg, nu, im, x = const_g2(key, e, Pmax)
        out = report(key, ps, lg, nu, e['c'], extra='(max|Im G| %.1e, %.0fs) %s' % (im, time.time() - t0, x if which != 'g2' else [b for b in x if b[0] != 2] if x else ''))
        res[key] = {str(k): v for k, v in out.items()}
    json.dump(res, open(D + '/s6_%s_%s.json' % (which, keys[0] if len(keys) == 1 else 'batch%d' % len(keys)), 'w'))
