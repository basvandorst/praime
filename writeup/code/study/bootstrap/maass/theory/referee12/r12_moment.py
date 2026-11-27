# Referee code, stage 3: the note's cosine-truncation lower bound, recomputed from my own moments,
# compared with the exact Delta_B; Lemma 5 probes; sign of A'_B(s) below s = 4; sign of g_B.
import mpmath as mp, json, sys
mp.mp.dps = 60
HERE = '/home/bas/praime/writeup/code/study/bootstrap/maass/theory/'
sys.path.insert(0, HERE + 'referee12')
from r12_main import G, loadB, setup_P, Pval, Pmin, xp, PR

def LB_Delta(Ms, k, mmax):
    # 2 M0 [psi(6) - psi(k/2)] + 2 sum_{j=1}^{2m+1} (-1)^j M_{2j} [zeta(2j+1,k/2) - zeta(2j+1,6)]
    a = mp.mpf(k)/2
    base = 2*Ms[0]*(mp.digamma(6) - mp.digamma(a))
    best, acc = None, base
    for j in range(1, 2*mmax+2):
        acc += 2*(-1)**j*Ms[j]*(mp.zeta(2*j+1, a) - mp.zeta(2*j+1, 6))
        if j % 2 == 1:
            best = acc if best is None else max(best, acc)
    return best

def LB_Aprime(Ms, s, mmax):
    a = mp.mpf(s)/2; best, acc = None, 0
    for j in range(0, 2*mmax+2):
        acc += (-1)**j*(2*j+1)*mp.zeta(2*j+2, a)*Ms[j]
        if j % 2 == 1:
            best = acc if best is None else max(best, acc)
    return best

EPS_NOTE = {0: '4.086', 8: '8.785e-3', 18: '3.694e-8'}
out = {}
for deg in (0, 8, 18):
    b = loadB(deg); Gb = G(24, b); g0 = Gb.g(0)
    Ms = [Gb.M(j) for j in range(0, 40)]
    P = setup_P(Gb)
    D = mp.fsum(Pval(P[p], xp(p)) - Pmin(P[p])[0] for p in PR)
    r = {}
    r['M0/g0 (should be 1)'] = mp.nstr(Ms[0]/g0, 20)
    r["LB A'(12)/g0"] = mp.nstr(LB_Aprime(Ms, 12, 19)/g0, 8)
    r["exact A'(12)/g0"] = mp.nstr(Gb.Aprime(12)/g0, 8)
    k = 12 - mp.mpf(EPS_NOTE[deg])
    lb = LB_Delta(Ms, k, 19); ex = Gb.Delta(k)
    r['at note kappa: LB Delta/(2D)'] = mp.nstr(lb/(2*D), 10)
    r['at note kappa: exact Delta/(2D)'] = mp.nstr(ex/(2*D), 10)
    r['LB <= exact'] = bool(lb <= ex)
    # reproduce the note's eps by bisection on the LB
    lo, hi = mp.mpf(0), mp.mpf(12)  # eps
    for _ in range(80):
        mid = (lo+hi)/2
        if LB_Delta(Ms, 12-mid, 19) > 2*D: hi = mid
        else: lo = mid
    r['eps from LB (bisection)'] = mp.nstr(hi, 6)
    # A'_B(s) on (0, 4]: sign (Lemma 5 does not cover it)
    r["min A'(s)/g0 on s in {0.05,0.2,0.5,1,2,3,4}"] = mp.nstr(min(Gb.Aprime(s) for s in [0.05, 0.2, 0.5, 1, 2, 3, 4])/g0, 6)
    # sign of g_B on [0, 6]
    us = [mp.mpf(i)/50 for i in range(0, 301)]
    gv = [Gb.g(u) for u in us]
    r['g_B sign changes on [0,6] (grid 0.02)'] = sum(1 for i in range(300) if gv[i]*gv[i+1] < 0)
    r['min g_B/g0 on grid'] = mp.nstr(min(gv)/g0, 4)
    if deg == 8:
        v, xm = Pmin(P[2]); r['deg8 P_2 argmin, Defect_2/g0'] = [mp.nstr(xm, 6), mp.nstr((Pval(P[2], xp(2))-v)/g0, 4)]
    out[deg] = r
    print(deg, json.dumps(r)); sys.stdout.flush()

# Lemma 5 probes: inf_y Re psi'(x+iy) on a grid, for x around and below 2
for x in ['2', '1.5', '1.2', '1.0', '0.8', '0.6', '0.5', '0.3']:
    x = mp.mpf(x)
    ys = [mp.mpf(i)/10 for i in range(0, 3001)] + [mp.mpf(10)**e for e in range(3, 7)]
    vals = [(mp.re(mp.psi(1, x+1j*y)), y) for y in ys]
    v, y = min(vals)
    print('Lemma5 probe x=%s: min Re psi\' = %s at y = %s' % (mp.nstr(x, 3), mp.nstr(v, 5), mp.nstr(y, 5)))
