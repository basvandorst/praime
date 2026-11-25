"""T5g: genus-2 curves of GL(2) type found in the T4b search (mean |lambda_p|^2 ~ 2): conductors 529 = 23^2 and 961 = 31^2.
Check that L(C) = L(f) L(f^sigma) for the weight-2 newforms of level 23 / 31 with real-quadratic coefficients (Dirichlet
convolution of the two embeddings' coefficients against lfunan of the genus-2 L, n <= 300), then attach the factor
constants frozen in round 1 (../predictions.json: W23a, W23b, W31a, W31b; their own zeros were never computed).
Readings (registered): half phase theta_f = theta_C/2 on all zeros: -(eps_a c_a + eps_b c_b) T2/2pi; product phase: 0.
usage: b2b_t5g.py <outdir>"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from bz_lib import build, params
from ulib import to_c, theta
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
TS = json.load(open('testset.json')); OUT = sys.argv[1]
P1 = json.load(open('../predictions.json'))['Lfunctions']
T2 = TS.get('height', {}).get('T5g', 150)
for N, keys in ((23, ('W23a', 'W23b')), (31, ('W31a', 'W31b'))):
    rows = [r for r in TS['T4b']['all_valid'] if r[0] == N * N]
    assert rows, N
    f, h = rows[0][1], rows[0][2]
    L = build(pari, ['g2', f, h]); Q, k, mu = params(pari, L)
    an = [to_c(x).real for x in pari.lfunan(L, 300)]
    mf = pari.mfinit([N, 2], 0); Lv = pari.lfunmf(mf, pari.mfeigenbasis(mf)[0])
    A = [[to_c(x).real for x in pari.lfunan(Lv[e], 300)] for e in (0, 1)]
    conv = [sum(A[0][d - 1] * A[1][n // d - 1] for d in range(1, n + 1) if n % d == 0) for n in range(1, 301)]
    err = max(abs(conv[i] - an[i]) for i in range(300))
    eps = to_c(pari.lfunrootres(L)[2])
    ca, cb = P1[keys[0]]['c'], P1[keys[1]]['c']; ea, eb = complex(*P1[keys[0]]['eps']), complex(*P1[keys[1]]['eps'])
    half = -(ea * ca + eb * cb) * T2 / (2 * math.pi)
    r = dict(key='T5g_N%d' % (N * N), cls='T5g', spec=['g2', f, h], Q=Q, k=k, mu=mu, degree=4, eps=[eps.real, eps.imag], T2=T2,
             M=float(theta(T2, Q, mu)) / math.pi, c=0.0, c_method='0 (L(f) L(f^sigma), verified: max |a_n - (a*b)_n| = %.1e, n <= 300)' % err,
             factors=dict(keys=keys, c=[ca, cb], eps=[[ea.real, ea.imag], [eb.real, eb.imag]], Q=N, mu=[0.5, 1.5]),
             half_phase_pred=[half.real, half.imag], feq=float(pari.lfuncheckfeq(L)), mean_lam2=rows[0][3])
    json.dump(r, open(os.path.join(OUT, r['key'] + '.json'), 'w'), indent=1)
    print(r['key'], 'factorisation error', err, 'c_a, c_b', ca, cb, 'half-phase prediction %.2f' % half.real, 'eps', eps, 'feq', r['feq'])
