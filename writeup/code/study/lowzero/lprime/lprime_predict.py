"""lprime_predict.py: frozen predictions for log|L'(rho, chi_D)| at the zeros 100 <= gamma <= 600 (derivation.md).
Writes predictions_lprime.json with, for every discriminant of the sets B, R120, Dn, Fn, the per-L-function predicted values
of every statistic of lprime_test.py, under the hypotheses
  LDM   leading order (recipe / hybrid Euler-Hadamard) + the first-order local-density-modulation term (1/L),
  lead  leading order only (no 1/L term),
  none  no dependence on the conductor's primes (every prime treated as good),
and the exact recipe (ratios conjecture, Conrey-Snaith sect. 7.1 with zeta -> zeta_q) for sum |L'(rho)|^2.
Only the discriminants and the heights 100..600 are used (zero-density weights; no zero of any set is read).
usage: lprime_predict.py OUT.json        (self-contained: numpy, mpmath)
"""
import sys, json, math
import numpy as np, mpmath as mp
mp.mp.dps = 30
EG = float(mp.euler); LO, HI = 100.0, 600.0
PR = (2, 3, 5, 7, 11, 13)


def is_fund(D):
    if D % 4 == 1: return all(D % (p * p) for p in range(3, int(abs(D) ** .5) + 1, 2))
    if D % 4 == 0:
        m = D // 4
        if m % 4 not in (2, 3): return False
        return all(m % (p * p) for p in range(3, int(abs(m) ** .5) + 1, 2))
    return False


def primes_of(n):
    f, p = [], 2
    while p * p <= n:
        if n % p == 0:
            f.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: f.append(n)
    return f


SETS = {'B': lambda D: 400 < abs(D) <= 700, 'R120': lambda D: 100 <= abs(D) <= 400,
        'Dn': lambda D: 700 < abs(D) <= 1000, 'Fn': lambda D: 1000 < abs(D) <= 1300}


def set_discriminants(name):
    lo = {'B': 401, 'R120': 100, 'Dn': 701, 'Fn': 1001}[name]; hi = {'B': 700, 'R120': 400, 'Dn': 1000, 'Fn': 1300}[name]
    return [D for a in range(lo, hi + 1) for D in (-a, a) if D not in (0, 1) and is_fund(D)]


# ---------- CUE(N) ----------
def lnM(n, s):
    G = mp.barnesg
    return mp.log(G(n + 1)) + mp.log(G(n + 2 * s + 1)) + 2 * mp.log(G(1 + s)) - 2 * mp.log(G(n + s + 1)) - mp.log(G(1 + 2 * s))


def cue_tail_C(N):
    """P(|Lambda'(theta_n)| < y) ~ C y^3 in CUE(N): E|Lambda'|^{2k} = M_{N-1}(1+k)/N ~ (3C/2)/(k + 3/2)."""
    h = mp.mpf('1e-12'); k = -1.5 + h
    return float(2 * h * mp.exp(lnM(N - 1, 1 + k)) / N / 3)


# ---------- one Euler factor ----------
MPH = 1 << 14; PHI = 2 * np.pi * (np.arange(MPH) + 0.5) / MPH


def _ser_mul(a, b, n): return np.array([sum(a[i] * b[j - i] for i in range(j + 1)) for j in range(n + 1)])


def _ser_inv(a, n):
    b = np.zeros(n + 1); b[0] = 1 / a[0]
    for j in range(1, n + 1): b[j] = -sum(a[i] * b[j - i] for i in range(1, j + 1)) / a[0]
    return b


def _ser_log(a, n):
    d = np.array([(j + 1) * a[j + 1] for j in range(n)]); q = _ser_mul(d, _ser_inv(a, n), n - 1)
    return np.concatenate([[0.0], [q[j - 1] / j for j in range(1, n + 1)]])


def euler_cumulants(p, order=4):
    """(lead, ldm): cumulants kappa_1..order of A_p (uniform phase) and the LDM first-order term per unit (ln p)/L:
    K1(t) = -2 (1 + t/2)^2 E[e^{tA} R]/E[e^{tA}], R = Re z/(1-z), z = p^{-1/2} e^{i phi}."""
    z = p ** -0.5 * np.exp(1j * PHI); A = -np.log(np.abs(1 - z)); R = np.real(z / (1 - z)); n = order
    def mom(w):
        out, fj, fa = [], np.ones(MPH), 1.0
        for j in range(n + 1): out.append(float(np.mean(w * fj)) / fa); fj = fj * A; fa *= j + 1
        return np.array(out)
    E0 = mom(np.ones(MPH)); ER = mom(R); K0 = _ser_log(E0, n); G = _ser_mul(ER, _ser_inv(E0, n), n)
    pre = np.zeros(n + 1); pre[0], pre[1], pre[2] = 1, 1, 0.25; K1 = -2 * _ser_mul(pre, G, n)
    fact = np.array([math.factorial(j) for j in range(n + 1)], float)
    return (K0 * fact)[1:], (K1 * fact)[1:]


def G_of(p, t):
    z = p ** -0.5 * np.exp(1j * PHI); w = np.abs(1 - z) ** (-t); R = np.real(z / (1 - z)); return float(np.mean(w * R) / np.mean(w))


def hyp(k, p): return float(mp.hyp2f1(k, k, 1, mp.mpf(1) / p))


# ---------- exact recipe for k = 1 ----------
_SZ = [mp.mpf(1)] + [(-1) ** j * mp.stieltjes(j) / mp.factorial(j) for j in range(6)]


def laurent(P, n=4):
    """c_0..c_{n-1} of Z(s) = zeta_q(1+s)/c_q = 1/s + c_0 + c_1 s + ..., zeta_q(1+s) = zeta(1+s) prod_{p|q} (1 - p^{-1-s})."""
    ser = _SZ[:n + 2]   # s*zeta(1+s)
    for p in P:
        lp = mp.log(p)   # (1 - p^{-1} p^{-s}) / (1 - 1/p) as a series in s
        fac = [(1 - mp.mpf(1) / p) if j == 0 else -(mp.mpf(1) / p) * (-lp) ** j / mp.factorial(j) for j in range(n + 2)]
        fac = [f / (1 - mp.mpf(1) / p) for f in fac]
        ser = [sum(ser[i] * fac[j - i] for i in range(j + 1)) for j in range(n + 2)]
    return [ser[j + 1] for j in range(n)]


def m1_poly(L, P):
    """Density in t of sum_gamma |L'(rho)|^2 (recipe): (c_q/2pi)[L^4/12 + (2/3)c0 L^3 + (c0^2 + 2c1) L^2
    + (-2c0^3 + 10 c0 c1 - 2 c2) L + 2c0^4 - 12 c0^2 c1 + 14 c1^2 + 16 c0 c2 - 20 c3]  (Conrey-Snaith Thm 7.1 rewritten in
    the Laurent coefficients of Z; for q = 1, c_n = (-1)^n gamma_n/n!)."""
    c0, c1, c2, c3 = [float(x) for x in laurent(P)]; cq = float(np.prod([1 - 1 / p for p in P])) if P else 1.0
    return cq / (2 * math.pi) * (L ** 4 / 12 + 2 / 3 * c0 * L ** 3 + (c0 ** 2 + 2 * c1) * L ** 2 + (-2 * c0 ** 3 + 10 * c0 * c1 - 2 * c2) * L
                                 + 2 * c0 ** 4 - 12 * c0 ** 2 * c1 + 14 * c1 ** 2 + 16 * c0 * c2 - 20 * c3)


# ---------- heights ----------
GX, GW = np.polynomial.legendre.leggauss(48)
TN = LO + (HI - LO) * (GX + 1) / 2; TW = (HI - LO) / 2 * GW


def density(t, q, a):
    """N'(t) of the smooth count theta_chi(t)/pi: (1/2pi)[ln(q/pi) + Re psi((1/2 + a + i t)/2)] (a = 0 even, 1 odd)."""
    return np.array([(math.log(q / math.pi) + float(mp.re(mp.digamma((0.5 + a + 1j * x) / 2)))) / (2 * math.pi) for x in t])


def per_lfunction(D, consts):
    q = abs(D); a = 0 if D > 0 else 1; P = primes_of(q); Lt = np.log(q * TN / (2 * math.pi)); rho = density(TN, q, a)
    n = float((rho * TW).sum()); avg = lambda f: float((rho * TW * f).sum() / n)
    inv = avg(1 / Lt); lnL = avg(np.log(Lt))
    out = {'D': D, 'bad': [p for p in P], 'n_zeros_expected': n, 'mean_lnL': lnL, 'mean_invL': inv}
    # per-prime contributions (bad primes remove the good-prime contribution)
    lead = {n_: 0.0 for n_ in (1, 2, 3, 4)}; ldm = {n_: 0.0 for n_ in (1, 2, 3, 4)}
    for p in P:
        K0, K1 = consts['euler'](p)
        for j in range(4): lead[j + 1] -= K0[j]; ldm[j + 1] -= K1[j] * math.log(p) * inv
    for j in (1, 2, 3, 4):
        out['cum%d_lead' % j] = lead[j]; out['cum%d_ldm' % j] = lead[j] + ldm[j]; out['cum%d_none' % j] = 0.0
    # moments: log of the bad-prime factor on E|L'|^{2k}
    for k in (-1.0, -0.5, 0.5, 1.0):
        lg_lead = sum(-math.log(hyp(k, p)) for p in P)
        lg_ldm = sum(-math.log(hyp(k, p)) + math.log(1 + 2 * math.log(p) * (1 + k) ** 2 * consts['G'](p, 2 * k) * inv) for p in P)
        out['logmom%+.1f_lead' % k] = lg_lead; out['logmom%+.1f_ldm' % k] = lg_ldm; out['logmom%+.1f_none' % k] = 0.0
    # exact recipe k = 1 (absolute): integral of m1 over [100, 600]; and its no-conductor-prime alternative
    m1 = np.array([m1_poly(x, P) for x in Lt]); m1n = np.array([m1_poly(x, []) for x in Lt])
    m1lead = np.array([np.prod([1 - 1 / p for p in P]) / (2 * math.pi) * x ** 4 / 12 for x in Lt])
    out['S1_recipe'] = float((m1 * TW).sum()); out['S1_none'] = float((m1n * TW).sum()); out['S1_lead'] = float((m1lead * TW).sum())
    # k = -1 (Pearce-Crump): sum over zeros of |L'|^-2 ~ A_q/2pi * (HI - LO); A_q = (6/pi^2) prod p/(p+1)
    Aq = 6 / math.pi ** 2 * float(np.prod([p / (p + 1) for p in P])) if P else 6 / math.pi ** 2
    out['Sm1_conj'] = Aq / (2 * math.pi) * (HI - LO); out['Sm1_none'] = 6 / math.pi ** 2 / (2 * math.pi) * (HI - LO)
    # left tail constant C_q(L) = a_{-3/2}(q) C_CUE(N = L) (leading order) and its LDM factor; tail term for truncation
    out['tail_logfactor_lead'] = sum(math.log(1 / hyp(-1.5, p)) for p in P)
    out['tail_logfactor_ldm'] = sum(math.log(1 / hyp(-1.5, p)) + math.log(1 + 2 * math.log(p) * 0.25 * consts['G'](p, -3.0) * inv) for p in P)
    aq = consts['a32'] * math.exp(out['tail_logfactor_lead']); out['tail_C_weighted'] = float((rho * TW * np.array([aq * consts['CN'](x) for x in Lt])).sum())
    return out


def main(path):
    eul = {p: euler_cumulants(p) for p in range(2, 1400) if all(p % r for r in range(2, int(p ** .5) + 1))}
    Gc = {}
    def Gf(p, t):
        if (p, t) not in Gc: Gc[(p, t)] = G_of(p, t)
        return Gc[(p, t)]
    # a_{-3/2}(1) = prod_p (1-1/p)^{9/4} 2F1(-3/2,-3/2;1;1/p) (converges like 1 + O(p^-2))
    a32 = 1.0
    for p in eul: a32 *= (1 - 1 / p) ** 2.25 * hyp(-1.5, p)
    Ns = np.linspace(6, 14, 17); Cs = [cue_tail_C(x) for x in Ns]
    consts = {'euler': lambda p: eul[p], 'G': Gf, 'a32': a32, 'CN': lambda x: float(np.interp(x, Ns, Cs))}
    # global constants
    cA = sum(0.5 * (float(mp.polylog(2, mp.mpf(1) / p)) + math.log(1 - 1 / p)) for p in eul) + sum(
        0.5 * (float(mp.polylog(2, mp.mpf(1) / p)) + math.log(1 - 1 / p)) for p in range(1400, 200000) if all(p % r for r in range(2, int(p ** .5) + 1)))
    C1 = sum(math.log(p) * sum((1 - 1 / k) / p ** k for k in range(2, 60)) for p in eul)
    J = {'version': 'lprime-predictions-1', 'heights': [LO, HI], 'primes_design': list(PR),
         'constants': {'cue_var_const': 0.5 * (EG + 1) + 1 - math.pi ** 2 / 4, 'cA_var_arith': cA, 'C1_mean': C1,
                       'a_-3/2(1)': a32, 'cue_tail_C(N)': dict(zip([float(x) for x in Ns], Cs))},
         'per_prime': {str(p): {'kappa_lead': list(eul[p][0]), 'kappa_ldm_per_lnp_over_L': list(eul[p][1]),
                                'r_k': {str(k): 1 / hyp(k, p) for k in (-1.5, -1, -0.5, 0.5, 1)},
                                'ldm_k': {str(k): 2 * math.log(p) * (1 + k) ** 2 * Gf(p, 2 * k) for k in (-1.5, -1, -0.5, 0.5, 1)}}
                       for p in PR}}
    for name in ('B', 'R120', 'Dn', 'Fn'):
        J[name] = [per_lfunction(D, consts) for D in set_discriminants(name)]
        print(name, len(J[name]), flush=True)
    json.dump(J, open(path, 'w'), indent=0)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--check':
        for P, L in (([], 10.0), ([2], 10.0), ([3], 10.0), ([], 5.0)):
            print('m1_poly', P, L, m1_poly(L, P), 'ratio to q=1' if P else '', m1_poly(L, P) / m1_poly(L, []) if P else '')
        print('cue tail C(N) for N = 8, 10, 12:', [cue_tail_C(N) for N in (8, 10, 12)])
        sys.exit()
    main(sys.argv[1])
