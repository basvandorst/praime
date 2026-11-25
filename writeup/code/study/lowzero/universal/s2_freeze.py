"""freeze the predictions of the universal one-point law (no zero is read by this script).
For each L in the catalogue: window (0, T2]; predicted sharp sum -eps c T2/2pi, Hann (sin^2 taper on (0,T2]) half of it;
lines r = 2, 1/2, 3, 1/3: -eps c S_p(v)/F_p(1/p) * T2/2pi; heuristic noise sd = sqrt(kappa M), kappa = 0.25 (degree-2 data
of round 148: rms 13.9 with M ~ 840 zeros), M = theta(T2)/pi; products: 0 and the cross-structure readings.
Writes predictions.json and prints its sha256."""
import sys, os, json, math, cmath, hashlib, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ulib import CATALOGUE, local_data, theta
KAPPA = 0.25
pred = dict(created_utc=datetime.datetime.utcnow().isoformat(timespec='seconds'),
            law='sum_{0<gamma<=T2} e^{2 i theta_L(gamma)} = -eps_L c_L T2/2pi + o(T), c_L = B_L(1)/L(1, Ad L); c_L = 0 for non-primitive L',
            conventions='theta_L(t) = (t/2) ln(Q/pi^d) + sum_j Im lnGamma((1/2+mu_j+it)/2); eps from PARI lfunrootres (Lambda(s) = eps conj-Lambda(1-s)); '
                        'Hann weight sin^2(pi gamma/T2) (prediction = half of the sharp one); lines: sum e^{2 i theta(gamma)} r^{-i gamma}',
            noise='heuristic sd = sqrt(0.25 M) (sharp), sqrt(0.25*3M/8) (Hann), M = theta(T2)/pi; the decision uses the 6-block standard error '
                  '(blocks of equal t-length), reported together with the heuristic',
            decisions={
                'per L (c != 0)': 'slope b = Re(G/pred) (sharp); PASS if |b - 1| <= 3 se_block; also report Hann slope and Im(G/pred) (expected 0)',
                'per family': 'pooled slope sum Re(G conj pred)/sum |pred|^2, se from the per-L block variances; PASS if within 3 se of 1',
                'phase (non-self-dual)': 'arg(G) against arg(-eps c); alternatives: -c (no root number), -conj(eps) c; report the angle in degrees with error se/|pred|',
                'products': 'effective c_hat = -G 2pi/(eps T2); PASS if |c_hat| <= 3 se; cross-structure readings below',
                'lines': 'pooled over all primitive L: slope of Re(G_r conj pred_r) on |pred_r|^2, se from the heuristic noise per line (sd as for r = 1)',
                'alternatives': 'H0: no main term (all slopes 0); H_noroot: -c T/2pi (without eps); H_conj: -conj(eps) c T/2pi; for products H_naive: c_direct(10^5) instead of 0'},
            Lfunctions={})
for key, spec, T2, fam in CATALOGUE:
    fn = 'constants_%s.json' % key
    if not os.path.exists(fn):
        print('missing', key); continue
    d = json.load(open(fn))
    eps = complex(*d['eps']); c = float(d['c']); Q = d['Q']; mu = d['mu']
    M = float(theta(T2, Q, mu)) / math.pi
    unit = T2 / (2 * math.pi)
    G = -eps * c * unit
    lines = {}
    for p in (2, 3):
        P = [complex(a, b) for a, b in d['poly_2_3'][str(p)]]
        F, S, lam = local_data(P, p)
        for v in (1, -1):
            amp = -eps * c * S[v] / F if c != 0 else 0j
            lines['%s%d' % ('' if v > 0 else '1/', p)] = dict(amp_per_unit=[amp.real, amp.imag], sharp=[(amp * unit).real, (amp * unit).imag],
                                                              rho_over_sqrt_ab=[(S[v] / F).real, (S[v] / F).imag])
    e = dict(family=fam, spec=d['spec'], Q=Q, degree=len(mu), mu=mu, eps=[eps.real, eps.imag], c=c, c_method=d['c_method'],
             c_direct_smooth=d['c_direct_smooth'][0], c_direct_end_1e5=d['c_direct_end'][0], T1=0.0, T2=T2, M_expected=M,
             pred_sharp=[G.real, G.imag], pred_hann=[G.real / 2, G.imag / 2],
             sd_heur_sharp=math.sqrt(KAPPA * M), sd_heur_hann=math.sqrt(KAPPA * 3 * M / 8), lines=lines,
             alt_noroot=[-c * unit, 0.0], alt_conj=[(-eps.conjugate() * c * unit).real, (-eps.conjugate() * c * unit).imag],
             mean_abs_lambda_p_sq=d['mean_abs_lambda_p_sq'])
    if fam.startswith('product'):
        e['alt_naive_direct'] = [(-eps * d['c_direct_end'][0] * unit).real, 0.0]
    pred['Lfunctions'][key] = e
# cross-structure readings for the products (zero sets split by the factor whose zero it is)
cr = {}
for Z, R in (('Z23', 'R23'), ('Z31', 'R31')):
    if Z in pred['Lfunctions'] and R in pred['Lfunctions']:
        T2 = pred['Lfunctions'][Z]['T2']; unit = T2 / (2 * math.pi); cR = pred['Lfunctions'][R]['c']
        cr[Z] = {'all zeros, phase theta_K = theta_zeta + theta_rho': 0.0,
                 'zeta zeros, phase theta_zeta (degree-1 theorem)': -unit,
                 'rho zeros, phase theta_rho (degree-2 law)': -cR * unit,
                 'zeta zeros, phase theta_K': 0.0, 'rho zeros, phase theta_K': 0.0,
                 'all zeros, phase theta_zeta': -unit, 'all zeros, phase theta_rho': -cR * unit}
if 'J13' in pred['Lfunctions'] and 'N13' in pred['Lfunctions']:
    T2 = pred['Lfunctions']['J13']['T2']; unit = T2 / (2 * math.pi)
    n = pred['Lfunctions']['N13']; eps = complex(*n['eps']); cN = n['c']
    cr['J13'] = {'all zeros, phase theta_J = 2 theta_f': 0.0,
                 'all zeros, half phase theta_f (union of the f and f-bar laws)': -(eps + eps.conjugate()).real * cN * unit,
                 'f zeros, phase theta_f': [(-eps * cN * unit).real, (-eps * cN * unit).imag],
                 'f-bar zeros, phase theta_f': [(-eps.conjugate() * cN * unit).real, (-eps.conjugate() * cN * unit).imag]}
pred['cross_structure'] = cr
s = json.dumps(pred, indent=1, sort_keys=True)
open('predictions.json', 'w').write(s)
print('L-functions frozen:', len(pred['Lfunctions']))
for k, e in pred['Lfunctions'].items():
    print('%-5s %-32s d=%d Q=%-6g eps=%+.4f%+.4fi c=%+.5f T2=%4d M=%6.0f pred=%+9.2f%+9.2fi sd~%.1f  lines 2:%+.2f 1/2:%+.2f 3:%+.2f 1/3:%+.2f' % (
        k, e['family'][:32], e['degree'], e['Q'], e['eps'][0], e['eps'][1], e['c'], e['T2'], e['M_expected'], e['pred_sharp'][0], e['pred_sharp'][1],
        e['sd_heur_sharp'], e['lines']['2']['sharp'][0], e['lines']['1/2']['sharp'][0], e['lines']['3']['sharp'][0], e['lines']['1/3']['sharp'][0]))
print(json.dumps(cr, indent=1))
print('sha256(predictions.json) =', hashlib.sha256(s.encode()).hexdigest())
