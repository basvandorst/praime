"""dry-run set: the 23 zero sets of round 1 (+ zeta) under the blind key names, with the registered completeness rule
applied to them (m0 by lfunorderzero, eps by lfunrootres); the dry test set for b2_constants; the T5g analogue X1(13)
(L(f)L(f-bar), factor constants from ../predictions.json). Two proxy pairs of round 144 (not test pairs)."""
import sys, os, json, math, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from bz_lib import build, params, check
from ulib import to_c, theta
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
TS = json.load(open('testset.json'))
dryts = dict(seed=TS['seed'], height=dict(T1=1000, T2=600, T3=200, T4a=70, T4b=200, T5=300, T5g=200, T5pairs=400),
             T1=dict(test=[[7, 2, 4], [5, 2, 3]]), T2=dict(test=[[k, 1, 0] for k in (12, 16, 18, 20, 22, 26)]),
             T3=dict(test=[[11, [0, -1, 1, -10, -20]], [14, [1, 0, 1, 4, -6]], [15, [1, 1, 1, -10, -10]], [17, [1, -1, 1, -1, -14]],
                           [19, [0, 1, 1, -9, -15]], [37, [0, 0, 1, -1, 0]]]),
             T4a=dict(test=[]), T4b=dict(test=[[249, 'x^2+x', 'x^3+1'], [277, 'x^5-2*x^3+2*x^2-x', '1']]),
             T5=dict(test=[[-23, 'x^3-x-1'], [-31, 'x^3+x-1']]), T5pairs=dict(test=TS['T5pairs']['proxy']))
json.dump(dryts, open('dry/testset_dry.json', 'w'), indent=1)
MAP = dict(X7=('T1_7_2', ['char', 7, 2], 0), X7d=('T1_7_4', ['char', 7, 4], 0), X5=('T1_5_2', ['char', 5, 2], 0), X5d=('T1_5_3', ['char', 5, 3], 0),
           S11a=('T3_N11_0_m1_1_m10_m20', ['sym', [0, -1, 1, -10, -20], 2], 0), S14a=('T3_N14_1_0_1_4_m6', ['sym', [1, 0, 1, 4, -6], 2], 0),
           S15a=('T3_N15_1_1_1_m10_m10', ['sym', [1, 1, 1, -10, -10], 2], 0), S17a=('T3_N17_1_m1_1_m1_m14', ['sym', [1, -1, 1, -1, -14], 2], 0),
           S19a=('T3_N19_0_1_1_m9_m15', ['sym', [0, 1, 1, -9, -15], 2], 0), S37a=('T3_N37_0_0_1_m1_0', ['sym', [0, 0, 1, -1, 0], 2], 0),
           G249=('T4b_N249_0', ['g2', 'x^2+x', 'x^3+1'], 0), G277=('T4b_N277_1', ['g2', 'x^5-2*x^3+2*x^2-x', '1'], 0),
           Z23=('T5_d23', ['nf', 'x^3-x-1'], 1), Z31=('T5_d31', ['nf', 'x^3+x-1'], 1), J13=('T5g_N169', ['g2', 'x^5+x^4', 'x^3+x+1'], 0),
           Zeta=('Zeta', ['zeta'], 1))
for k in (12, 16, 18, 20, 22, 26):
    MAP['D12' if k == 12 else 'F%d' % k] = ('T2_k%d_N1_0' % k, ['mf', k, 1, 0], 0)
for old, (new, spec, poles) in MAP.items():
    z = np.load('../zeros/zeros_%s.npz' % old); g = z['g']
    L = build(pari, spec); Q, k, mu = params(pari, L)
    eps = to_c(pari.lfunrootres(L)[2]); phi = 0.0 if abs(eps.imag) < 1e-9 else math.atan2(eps.imag, eps.real)
    m0 = int(pari.lfunorderzero(L)) if spec[0] != 'zeta' else 0
    cexp = 0.5 + poles - m0 / 2 - phi / (2 * math.pi)
    drift, rmean, ok = check(g, Q, mu, cexp)
    info = dict(spec=spec, Q=Q, mu=mu, m0=m0, c_exp=cexp, n=len(g), drift=drift, r_mean=rmean, status='COMPLETE' if ok else 'REJECTED', source='round 1 zeros_%s.npz' % old)
    np.savez_compressed('dry/zeros/%s.npz' % new, g=g, info=json.dumps(info))
    print('%-5s -> %-26s m0 %d c_exp %+.3f r_mean %+.3f drift %+.3f %s' % (old, new, m0, cexp, rmean, drift, info['status']), flush=True)
P1 = json.load(open('../predictions.json'))['Lfunctions']
n = P1['N13']; e = complex(*n['eps']); c = n['c']; T2 = 200
L = build(pari, ['g2', 'x^5+x^4', 'x^3+x+1']); Q, k, mu = params(pari, L)
half = -(e * c + e.conjugate() * c) * T2 / (2 * math.pi)
json.dump(dict(key='T5g_N169', cls='T5g', spec=['g2', 'x^5+x^4', 'x^3+x+1'], Q=Q, k=k, mu=mu, degree=4, eps=[1.0, 0.0], T2=T2,
               M=float(theta(T2, Q, mu)) / math.pi, c=0.0, c_method='0 (L(f) L(f-bar))', factors=dict(keys=['N13', 'N13d'], c=[c, c], eps=[[e.real, e.imag], [e.real, -e.imag]]),
               half_phase_pred=[half.real, half.imag]), open('dry/constants/T5g_N169.json', 'w'), indent=1)
print('dry T5g half-phase prediction %.2f' % half.real)
