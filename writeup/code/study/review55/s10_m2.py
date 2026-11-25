"""M2 = sum over rho-zeros (<= 300) of e^{2i theta_zeta}: own values, own -psi_rho(T/2pi), and the density term at the stationary
point t = 2pi of theta_zeta: (1/pi) int theta_rho'(t) e^{2i theta_zeta(t)} dt (computed numerically on (0, 300])."""
import json, math, os
import numpy as np, cypari2
from scipy.special import loggamma
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'; D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json')); zz = np.load(D + '/zeta_zeros_310.npy')
pari = cypari2.Pari(); pari.allocatemem(200 * 10 ** 6)
def theta(t, Q, mu):
    t = np.asarray(t, float); s = t / 2 * math.log(Q / math.pi ** len(mu))
    for m in mu: s = s + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return s
T2 = 300.0; x = T2 / (2 * math.pi)
tt = np.linspace(1e-6, T2, 600001); h = tt[1] - tt[0]
rows = []
for key, e in P['products'].items():
    if e['cls'] != 'T5': continue
    r = e['rho']; g = np.load(B + '/zeros_blind/' + key + '.npz')['g']
    isz = np.array([np.min(np.abs(zz - v)) < 1e-6 for v in g]); gr = g[~isz]; gr = gr[(gr > 0) & (gr <= T2)]
    M2 = complex(np.exp(2j * theta(gr, 1.0, [0.0])).sum())
    # psi_rho(x) from splitting: Lambda_rho(p^k) = ln p (sum_{f_i | k} f_i - 1)
    nf = pari.nfinit(pari(e['spec'][1])); psi = 0.0
    for p in [int(v) for v in pari('primes([2,%d])' % int(x))]:
        fs = [int(pr[3]) for pr in pari.idealprimedec(nf, p)]; k = 1; pk = p
        while pk <= x:
            psi += math.log(p) * (sum(f for f in fs if k % f == 0) - 1); k += 1; pk *= p
    # density term: (1/pi) int theta_rho' e^{2i theta_zeta} dt
    th_r = theta(tt, r['Q'], r['mu']); dth = np.gradient(th_r, h)
    dens = complex(np.sum(dth * np.exp(2j * theta(tt, 1.0, [0.0]))) * h / math.pi)
    rows.append((key, e['d'], M2, -psi, e['readings']['M2']['pred'], dens))
m = np.array([r[2].real for r in rows]); mech = np.array([r[3] for r in rows]); reg = np.array([r[4] for r in rows]); dens = np.array([r[5].real for r in rows])
print('own -psi_rho(T/2pi) vs registered mechanism: max |diff| %.2e' % np.max(np.abs(mech - reg)))
print('own M2 vs results_blind: max |diff| %.2e' % max(abs(r[2] - complex(*json.load(open(B + '/results_blind.json'))['products'][r[0]]['M2'])) for r in rows))
for name, pred in (('registered -psi', mech), ('-psi + density term', mech + dens)):
    b = np.sum(m * pred) / np.sum(pred ** 2); res = m - pred
    X = np.stack([np.ones(len(m)), pred], 1); bv, rss, *_ = np.linalg.lstsq(X, m, rcond=None); s2 = rss[0] / (len(m) - 2); cov = s2 * np.linalg.inv(X.T @ X)
    n = len(m); loo = np.array([(np.sum(m * pred) - m[i] * pred[i]) / (np.sum(pred ** 2) - pred[i] ** 2) for i in range(n)])
    print('%-22s slope %.3f (jack se %.3f) | residual m - pred: mean %+.2f rms %.2f | OLS a = %+.2f ± %.2f, b = %.3f ± %.3f' % (
        name, b, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2)), res.mean(), math.sqrt(np.mean(res ** 2)), bv[0], math.sqrt(cov[0, 0]), bv[1], math.sqrt(cov[1, 1])))
print('density term: mean %.2f, range %.2f..%.2f; mean ln|d| %.2f; Im parts mean %.2f' % (dens.mean(), dens.min(), dens.max(), np.mean([math.log(-r[1]) for r in rows]), np.mean([r[5].imag for r in rows])))
