"""(a) registered J13 reading that needs only the J13 zeros: sum over all zeros of e^{2i theta_f}, theta_f = theta_J/2 (frozen: -2 Re(eps_f) c_f T2/2pi);
(b) POST HOC mechanism for products (derived after the Z23 result; not registered): a degree-1 factor's chirp e^{2i theta_zeta}
    sweeps the frequency ln(t/2pi) and meets every integer line with aligned phase (stationary value e(-k) = 1):
    sum_{rho zeros<=T} e^{2i theta_zeta}            ~ -psi_rho(T/2pi)            (Landau lines of the rho-zero density; ADR's mechanism)
    sum_{rho zeros<=T} e^{2i theta_zeta + 2i theta_rho} ~ -c_rho sum_{k<=T/2pi} rho(1/k)  (lines of the one-point spectrum)
    the degree-2 chirp e^{2i theta_rho} on zeta zeros meets lines with phases e(-2 sqrt(k/Q)): no alignment, ~0."""
import sys, os, json, math, cmath, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from ulib import *
raw = open('predictions.json').read(); assert hashlib.sha256(raw.encode()).hexdigest() == sys.argv[1]
P = json.loads(raw)
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
out = {}
# (a)
e = P['Lfunctions']['J13']; T2 = e['T2']
g = np.load('zeros/zeros_J13.npz')['g']; g = g[(g > 0) & (g <= T2)]
thJ = theta(g, e['Q'], e['mu'])
half = complex(np.exp(1j * thJ).sum())
pred_half = P['cross_structure']['J13']['all zeros, half phase theta_f (union of the f and f-bar laws)']
blocks = [complex(np.exp(1j * thJ[(g > a) & (g <= a + T2 / 6)]).sum()) for a in np.linspace(0, T2 - T2 / 6, 6)]
pb = pred_half / 6
s_k = np.array([(b / pb).real for b in blocks])
out['J13_half_phase'] = dict(pred=pred_half, measured=[half.real, half.imag], slope=(half / pred_half).real,
                             se_block=float(s_k.std(ddof=1) / math.sqrt(6)))
print('(a) J13, all zeros, half phase theta_f: predicted %.1f, measured %.1f%+.1fi, slope %.3f ± %.3f (block)' % (
    pred_half, half.real, half.imag, (half / pred_half).real, s_k.std(ddof=1) / math.sqrt(6)))
# (b)
R = json.load(open('results.json'))
for Z, Rk in (('Z23', 'R23'), ('Z31', 'R31')):
    spec = [c for c in CATALOGUE if c[0] == Rk][0][1]
    L = build(pari, spec); Q, k, mu = params(pari, L)
    T2 = P['Lfunctions'][Z]['T2']; K = int(T2 / (2 * math.pi))
    an = pari.lfunan(L, K)
    # psi_rho(x) = sum_{p^v <= x} (alpha^v + beta^v) ln p from local polynomials
    psi = 0.0; rho = {}
    for p in primes_upto(K):
        p = int(p); Pp = local_poly(pari, L, p, k)
        roots = np.roots(Pp[::-1]) if len(Pp) > 1 else []
        alphas = [1 / r for r in roots]
        v = 1
        while p ** v <= T2 / (2 * math.pi):
            psi += sum(a ** v for a in alphas).real * math.log(p); v += 1
        lam = series_inv(Pp, 80)
        F = sum(np.conj(lam[j]) * Pp[j] * p ** (-j) for j in range(len(Pp)))
        for vv in range(1, 8):
            if p ** vv > K: break
            s = sum(np.conj(lam[h - vv]) * Pp[h] * p ** (-(2 * h - vv) / 2) for h in range(vv, len(Pp)))
            rho[(p, vv)] = p ** (vv / 2) * s / F
    M = 0
    for kk in range(1, K + 1):
        f = pari.factor(kk); val = 1
        for i in range(len(f[0])): val *= rho[(int(f[0][i]), int(f[1][i]))]
        M += val
    cR = P['Lfunctions'][Rk]['c']
    cr = R['cross'][Z]
    m_zeta_phase = complex(*cr['all zeros, phase theta_zeta']) - complex(*cr['zeta zeros, phase theta_zeta (degree-1 theorem)'])
    m_prod = complex(*cr['rho zeros, phase theta_K'])
    m_all = complex(*cr['all zeros, phase theta_K = theta_zeta + theta_rho'])
    out[Z] = dict(K=K, psi_rho=psi, sum_rho_1_over_k=M.real, mech_rho_zeros_zeta_phase=-psi, meas_rho_zeros_zeta_phase=[m_zeta_phase.real, m_zeta_phase.imag],
                  mech_rho_zeros_product_phase=-cR * M.real, meas_rho_zeros_product_phase=[m_prod.real, m_prod.imag],
                  meas_all_product_phase=[m_all.real, m_all.imag])
    print('(b) %s: x = T2/2pi = %.1f; rho zeros, zeta phase: mechanism -psi_rho(x) = %+.1f, measured %+.1f%+.1fi | rho zeros, product phase: mechanism -c sum rho(1/k) = %+.1f, measured %+.1f%+.1fi | all zeros product phase measured %+.1f%+.1fi (registered 0)' % (
        Z, T2 / (2 * math.pi), -psi, m_zeta_phase.real, m_zeta_phase.imag, -cR * M.real, m_prod.real, m_prod.imag, m_all.real, m_all.imag))
# phase angles of the degree-1 non-self-dual sums
for key in ('X7', 'X7d', 'X5', 'X5d'):
    r = R['per_L'][key]; e = P['Lfunctions'][key]
    G = complex(*r['G']); pr = complex(*e['pred_sharp'])
    out[key + '_angle'] = dict(arg_measured_deg=math.degrees(cmath.phase(G)), arg_pred_deg=math.degrees(cmath.phase(pr)),
                               arg_noroot_deg=180.0, arg_conj_deg=math.degrees(cmath.phase(complex(*e['alt_conj']))),
                               angle_err_deg=math.degrees(e['sd_heur_sharp'] / abs(pr)))
    print('    %s: arg G = %.1f deg; law (-eps c) %.1f; no root number 180.0; conj(eps) %.1f; heuristic error %.1f deg' % (
        key, math.degrees(cmath.phase(G)), math.degrees(cmath.phase(pr)), math.degrees(cmath.phase(complex(*e['alt_conj']))), math.degrees(e['sd_heur_sharp'] / abs(pr))))
json.dump(out, open('s5_extra.json', 'w'), indent=1, default=lambda z: [z.real, z.imag] if isinstance(z, complex) else str(z))
