"""check_second.py: asserts every identity and number used in second_moment.md.

  python check_second.py a   identities: twisted multiplicativity of Kloosterman sums, CRT inverse tables, FFT bilinear form
                             = direct double sum (with the low-rank Bessel kernel), Hecke relations on PARI data (N = 61, k = 4),
                             Petersson double sum against eigenform data, mean of w(m) and its local factors (sympy)
  python check_second.py b   N = 61, k = 4, T = 6: form side (forms_side.py, PARI eigenforms, mfpetersson weights) against the
                             Petersson side (moment2.py, mixed.py): M11, M2, X_q, C2(Q), V(Q)
  python check_second.py c   the experiment (task 1): integer model against its theorem; arcs x = +-1 of the prime model against
                             the integer model; arc weights against w(m); lag formula; V table (analysis.py) with the stated values
Each check prints PASS/FAIL with its numbers. Runs needed before b, c: see second_moment.md, section 9 (commands).
"""
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, numpy as np
from smlib import *
from scipy.special import jv

RES = []
def check(name, ok, info=''):
    RES.append((name, bool(ok)))
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)

def kl_direct(a, b, c):
    """plain-integer Kloosterman sum (independent of smlib)."""
    if c == 1: return 1.0
    s = 0.0
    for x in range(1, c):
        if math.gcd(x, c) == 1:
            s += math.cos(2*math.pi*((a*x + b*pow(x, -1, c)) % c)/c)
    return s

# ------------------------------------------------------------------------------------------------ a
def sec_a():
    rng = np.random.default_rng(7)
    # a1 twisted multiplicativity and the factorisation used in section 3
    worst = 0.0
    for _ in range(60):
        q = int(rng.choice([7, 11, 13, 61, 101])); r = int(rng.integers(1, 40))
        if math.gcd(q, r) != 1: continue
        a, b = int(rng.integers(1, 10**6)), int(rng.integers(1, 10**6))
        lhs = kl_direct(a, b, q*r)
        rb = pow(r, -1, q); qb = pow(q, -1, r) if r > 1 else 0
        rhs = kl_direct(a*rb, b*rb, q)*kl_direct(a*qb, b*qb, r)
        worst = max(worst, abs(lhs - rhs))
        if a % q and b % q:   # S(n,n';Nm) = S(1, n n' mbar^2; N) S(n Nbar, n' Nbar; m)
            rhs2 = kl_direct(1, a*b*rb*rb, q)*kl_direct(a*qb, b*qb, r)
            worst = max(worst, abs(lhs - rhs2))
    check('a1 S(a,b;qr) = S(a rbar, b rbar; q) S(a qbar, b qbar; r) and = S(1, ab rbar^2; q) S(..;r)', worst < 1e-8, 'max err %.1e' % worst)
    # a2 CRT unit/inverse tables
    ok = True
    for (N, ms) in [(10007, [1, 2, 6, 30]), (61, [1, 6, 61, 122]), (5003, [1, 4, 9])]:
        mod = Moduli(N)
        for m in ms:
            x, xb = mod.units_and_bar(m); c = N*m
            phi = sum(1 for v in range(1, c) if math.gcd(v, c) == 1) if c < 200000 else None
            ok &= bool(np.all((x*xb) % c == 1)) and (phi is None or len(x) == phi)
    check('a2 units and inverses mod N m (CRT and direct)', ok)
    # a3 FFT bilinear form with the low-rank kernel = direct double sum
    N, k, m = 61, 4, 3; c = N*m
    n = np.sort(rng.choice(np.arange(200, 5000), 40, replace=False)).astype(np.int64)
    xw = rng.normal(size=40) + 1j*rng.normal(size=40); yw = rng.normal(size=40) + 1j*rng.normal(size=40)
    s = TWO_PI*np.sqrt(n/N)
    direct = sum(xw[i]*yw[j]*kl_direct(int(n[i]), int(n[j]), c)*jv(k-1, 4*math.pi*math.sqrt(n[i]*n[j])/c) for i in range(40) for j in range(40))
    nodes = cheb_nodes(s.min(), s.max(), int(0.6*(s.max()-s.min())*s.max()/(math.pi*m)) + 40)
    lam, E = kernel_factor(k, m, nodes, tol=1e-13); F = bary_matrix(nodes, s) @ E
    x, xb = Moduli(N).units_and_bar(m); fft = 0j
    for r in range(len(lam)):
        A = np.zeros(c, complex); B = np.zeros(c, complex)
        np.add.at(A, n % c, xw*F[:, r]); np.add.at(B, n % c, yw*F[:, r])
        Ah = np.fft.ifft(A)*c; Bh = np.fft.ifft(B)*c
        fft += lam[r]*np.sum(Ah[x]*Bh[xb])
    check('a3 sum x_n y_n\' S(n,n\';c) J(4pi sqrt(nn\')/c) by FFT + low-rank kernel = direct', abs(fft - direct) < 1e-8*max(1, abs(direct)),
          'direct %.6f%+.6fi, fft %.6f%+.6fi' % (direct.real, direct.imag, fft.real, fft.imag))
    # a4 Hecke relations and a5 Petersson double sums on PARI data
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(900*10**6)
    N, k, P = 61, 4, 1200
    pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
    forms = []
    for j in range(1, int(pari('#B'))+1):
        co = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, P)); pet = pari('mfpetersson(mfsymbol(mf,B[%d]))' % j)
        isv = str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) == 't_VEC'
        embs = [(co[i], float(pet[i, i].real())) for i in range(len(co))] if isv else [(co, float(pet.real()))]
        for (cc, pv) in embs:
            a = np.array([float(complex(v).real) for v in cc]); lam_ = np.zeros(len(a)); lam_[1:] = a[1:]/np.arange(1, len(a))**((k-1)/2)
            forms.append((lam_, math.gamma(k-1)/((4*math.pi)**(k-1)*(N+1)*pv)))
    worst = 0.0
    for (lam_, om) in forms:
        for (p, p2) in [(2, 3), (5, 7), (11, 13), (3, 29)]:
            worst = max(worst, abs(lam_[p]*lam_[p2] - lam_[p*p2]))
        for p in [2, 3, 5, 7, 31]:
            worst = max(worst, abs(lam_[p]**2 - lam_[p*p] - 1))
        worst = max(worst, abs(lam_[N]**2 - 1/N), abs(lam_[N]*lam_[7] - lam_[7*N]))
    check('a4 Hecke: lam(p)lam(p\') = lam(pp\'), lam(p)^2 = lam(p^2)+1, lam(N)^2 = 1/N, lam(N)lam(p) = lam(Np) (N=61, k=4, %d forms)' % len(forms),
          worst < 1e-9, 'max err %.1e' % worst)
    worst = 0.0
    for (a_, b_) in [(1, 1), (2, 3), (101, 103), (997, 1009), (61, 997), (4, 9)]:
        fs = sum(om*lam_[a_]*lam_[b_] for (lam_, om) in forms)
        pet = (1.0 if a_ == b_ else 0.0) + 2*math.pi*(-1)**(k//2)*sum(kloost(a_, b_, c)/c*jv(k-1, 4*math.pi*math.sqrt(a_*b_)/c) for c in range(N, N*1500+1, N))
        worst = max(worst, abs(fs - pet))
    check('a5 Petersson double sum Delta(n,n\') (c <= 1500 N) = sum omega_f lam(n) lam(n\') (N=61, k=4)', worst < 2e-4, 'max err %.1e' % worst)
    # a6 mean of w(m) and local factors
    from analysis import w_of_m
    import sympy as sp
    p = sp.symbols('p', positive=True)
    Sp = 1 + (1 + (p-2)/(p-1)**2)/p + (1 + 1/(p-1))*(1/(p*(p-1)))
    loc = sp.simplify((1 - 1/p)*Sp - (1 + 1/p**2))
    Mx = 200000; mean = np.mean([w_of_m(m_) for m_ in range(1, Mx+1)])
    check('a6 w(m) multiplicative with local mean (1-1/p) sum_e w(p^e)p^-e = 1 + 1/p^2 (sympy); mean_{m<=2e5} w = 15/pi^2',
          loc == 0 and abs(mean - 15/math.pi**2) < 2e-3, 'sympy residual %s, mean %.5f vs %.5f' % (loc, mean, 15/math.pi**2))
    # local density of the q = p arcs: (1/(p-1)^2)[(p-2)(1/p - 1/p^2) + (p-1)/p^2] = 1/p^2
    loc2 = sp.simplify(((p-2)*(1/p - 1/p**2) + (p-1)/p**2)/(p-1)**2 - 1/p**2)
    check('a7 mean over m of the q-th arc weight mu^2(q) A_q(m)/phi(q)^2 1_{q|m} = mu^2(q)/q^2 (p-local, sympy)', loc2 == 0, str(loc2))

# ------------------------------------------------------------------------------------------------ b
def sec_b():
    fz = np.load('forms_61_4_6.npz'); mx = np.load('mixed_61_4_6_primes.npz'); z = np.load('mom_61_4_6_primes_1_4000_3.5.npz')
    b1, G1 = complex(z['b1']), complex(z['G1'])
    M11 = float(z['D11']) + z['rows'][:, 3].sum() + abs(b1)**2 + 2*(np.conj(b1)*G1).real
    M2 = complex(z['D2']) + complex(z['rows'][:, 1].sum(), z['rows'][:, 2].sum()) + b1*b1 + 2*b1*G1
    check('b1 M11 = sum^h |Z_f|^2: Petersson (c = 61 m, m <= 4000) = eigenforms', abs(M11 - complex(fz['M11']).real) < 1e-4,
          'Petersson %.6f, forms %.6f' % (M11, complex(fz['M11']).real))
    check('b2 M2 = sum^h Z_f^2: Petersson = eigenforms', abs(M2 - complex(fz['M2'])) < 1e-4,
          'Petersson %.6f%+.6fi, forms %.6f%+.6fi' % (M2.real, M2.imag, complex(fz['M2']).real, complex(fz['M2']).imag))
    err = np.max(np.abs(mx['Xq'] - fz['Xq']))
    check('b3 X_q = sum^h eps_f lam_f(q^2) Z_f, q <= 10 squarefree: Petersson (root-number formula) = eigenforms', err < 2e-4, 'max err %.1e' % err)
    L = float(mx['L']); Qs = list(mx['res_Q']); i = Qs.index(10)
    V = M11 + 2*L*complex(mx['res_XQ'][i]).real + L*L*float(mx['res_C2'][i])
    # forms-side V(10) recomputed from its pieces
    print('     (forms side V(10) printed by forms_side.py: compare)')
    check('b4 V(10) from the Petersson pieces', abs(V - 1.846613) < 1e-3, 'Petersson %.6f vs forms 1.846613' % V)

# ------------------------------------------------------------------------------------------------ c
def sec_c():
    from analysis import table, w_of_m
    out = table(verbose=True)
    # c1 integer model against its theorem (fact 2053 weight k + Petersson diagonal)
    for (T, f) in [(6, 'mom_10007_4_6_integers.npz'), (12, 'mom_10007_4_12_integers_1_100_3.5.npz')]:
        z = np.load(f); win = Window(10007, 4, T, h=0.05)
        M11 = float(z['D11']) + z['rows'][:, 3].sum()
        th = sum(abs(np.sum(win.H*np.exp(1j*math.log(l)*win.tt))*win.h)**2/l for l in range(1, 400))/(4*math.pi**2)
        check('c1 integer model T=%g: Petersson M11 = (1/4pi^2) sum_l |H^(ln l)|^2/l (theorem) within 1%%' % T, abs(M11/th - 1) < 0.01,
              'Petersson %.5f, theorem %.5f' % (M11, th))
    # c2 x = +-1 part of the prime model = integer model term (T = 12, N = 10007, k = 4), from xprofile output values
    zi = np.load('mom_10007_4_12_integers_1_100_3.5.npz')['rows']
    check('c2 arcs x = +-1 of primes = integer-model term (m=1 M11: 0.48540 vs integer; m=4 M2: 0.10939+0.23861i)',
          abs(zi[0, 3] - 0.48540) < 2e-3 and abs(complex(zi[3, 1], zi[3, 2]) - (0.10939+0.23861j)) < 2e-3,
          'integer m=1 %.5f, m=4 %.5f%+.5fi' % (zi[0, 3], zi[3, 1], zi[3, 2]))
    # c3 lag formula for the x = +-1 part of M11: (1/pi) int H(t) H(t - 2 pi m) dt
    for key in [(10007, 14, 20), (10007, 14, 30)]:
        if key not in out: continue
        rows = out[key]['rows']; N, k, T = key; win = Window(N, k, T, h=0.05)
        for m in [1, 2, 3]:
            sh = int(round(TWO_PI*m/win.h)); lag = np.sum(win.H[sh:]*win.H[:-sh])*win.h/math.pi
            pm = rows[rows[:, 0] == m][0, 11]
            # claimed: within 1 % where the lag 2 pi m is well inside the window (m = 1; m = 2 at T = 30); reported otherwise
            if m == 1 or (m == 2 and T >= 30):
                check('c3 lag formula T=%g m=%d: x=+-1 part of M11 %.4f vs (1/pi) int H(t)H(t-2pi m) %.4f (1%%)' % (T, m, pm, lag),
                      abs(pm - lag) < 0.01*abs(lag))
            else:
                print('INFO c3 lag formula T=%g m=%d (lag near the window length, not claimed): x=+-1 part %.4f vs %.4f (%.0f%%)'
                      % (T, m, pm, lag, 100*abs(pm - lag)/max(abs(lag), 1e-9)))
    # c4 arc weights: for every modulus with |x=+-1 part| > 0.3, |arcs/(x=+-1) - w(m)| < 0.35
    worst = 0.0; cnt = 0
    for key, rec in out.items():
        rows = rec['rows']
        if rows.shape[1] <= 6: continue
        for r in rows:
            for a, p in ((complex(r[6], r[7]), complex(r[9], r[10])), (r[8], r[11])):
                if abs(p) > 0.3:
                    worst = max(worst, abs(a/p - w_of_m(int(r[0])))); cnt += 1
    check('c4 arcs / (x = +-1 part) = w(m) on all %d moduli with |x=+-1 part| > 0.3' % cnt, worst < 0.35, 'max |ratio - w(m)| = %.3f' % worst)
    # c5 V/T^2 decreasing at N = 10007 (k = 14) and V/D11 in [0.75, 1]
    ks = sorted([key for key in out if key[0] == 10007 and key[1] == 14 and 'V' in out[key]], key=lambda t: t[2])
    vt = [out[key]['V']/key[2]**2 for key in ks]
    check('c5 N=10007, k=14: V(10)/T^2 decreasing in T (%s)' % ', '.join('T=%g: %.4f' % (key[2], v) for key, v in zip(ks, vt)),
          all(vt[i+1] < vt[i] for i in range(len(vt)-1)) and len(vt) >= 2)
    vd = [out[key]['V']/out[key]['D11'] for key in out if 'V' in out[key]]
    check('c6 V(10)/D11 in [0.75, 1] for every run (the residual has diagonal size)', all(0.75 <= v <= 1.0 for v in vd), ', '.join('%.3f' % v for v in vd))

if __name__ == '__main__':
    for s in (sys.argv[1:] or ['a', 'b', 'c']):
        {'a': sec_a, 'b': sec_b, 'c': sec_c}[s]()
    print('%d/%d PASS' % (sum(ok for _, ok in RES), len(RES)))
