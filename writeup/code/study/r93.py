"""Round 93 (preregister/round93.txt): the long-interval prime variance R(X, d) of rounds 91-92 from the zeros alone.
psi((1+d)x) - psi(x) - dx = -sum_rho c_rho x^rho, c_rho = ((1+d)^rho - 1)/rho (RH, explicit formula), so
R_zero(X, d) = (2/(d X^2)) Re sum_{rho,rho'} c_rho conj(c_rho') (X^(2+i(g-g')) - 1)/(2 + i(g-g')),
over the zeros up to height 2949 and their conjugates, plus the diagonal tail beyond. usage: r93.py"""
import os, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
Z = np.array([float(l) for l in open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'zeros_hp2419.txt')) if l.strip()])
T = float(Z[-1]); GAM = np.concatenate([Z, -Z]); RHO = 0.5 + 1j * GAM
G0 = 0.5772156649015329


def r_zero(X, d):
    c = ((1 + d) ** RHO - 1) / RHO
    lnX = math.log(X); diag = float(np.sum(np.abs(c) ** 2))
    tot = 0.0
    for i0 in range(0, len(GAM), 600):
        dg = GAM[i0:i0 + 600, None] - GAM[None, :]
        num = np.exp((2 + 1j * dg) * lnX) - 1.0
        tot += float(np.real(np.sum(c[i0:i0 + 600, None] * np.conj(c[None, :]) * num / (2 + 1j * dg))))
    tail = (2 + d) * (math.log(T / (2 * math.pi)) + 1) / (math.pi * T * d)
    R = 2 * tot / (d * X * X) + tail
    return {'R_zero': R, 'diagonal': diag / d + tail, 'offdiagonal': R - diag / d - tail, 'tail': tail}


if __name__ == '__main__':
    res = {'T': T, 'zeros': len(Z), 'rows': []}
    for d in (0.1, 0.01):
        R_hl = math.log(1 / d) - G0 - math.log(2 * math.pi) + (1 + d) * math.log(1 + d) / d
        for X in (1e8, 1e9, 1e10):
            data = json.load(open(os.path.join(OUT, 'r91_X%.0e.json' % X)))
            Rd = [r['R'] for r in data['rows'] if abs(r['delta'] - d) < 1e-12][0]
            z = r_zero(X, d)
            row = {'delta': d, 'X': X, 'R_primes': Rd, **z, 'R_HL_exact': R_hl, 'diff_zero_minus_primes': z['R_zero'] - Rd}
            res['rows'].append(row)
            print('d %.2f X %.0e  primes %.4f  zeros %.4f (diag %.4f + pairs %+.4f, tail %.4f)  diff %+.4f  HL %.4f' % (d, X, Rd, z['R_zero'], z['diagonal'], z['offdiagonal'], z['tail'], z['R_zero'] - Rd, R_hl), flush=True)
    json.dump(res, open(os.path.join(OUT, 'r93.json'), 'w'), indent=1)
