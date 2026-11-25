"""r3_const.py (review, own code): the constants c for A4, B6, C23, C47a, three ways:
 (1) ray form, conditionally convergent, natural order: c_P = (1/L(1,chi_D)) prod_{p<=P split} (1 - lam_{psi^2}(p)/(p-1))
 (2) the direct Euler product prod_{p<=P} c_p (c_p = 1 - (lam^2-1)/p split, 1 + 1/p inert, 1 ramified)
 (3) accelerated: (1/L(1,chi_D)) (1/L(1,psi^2)) prod_p E_p (absolutely convergent), L(1,psi^2) from
     A4: PARI lfun of my own coefficient vector (CM form, weight 9, level 4), B6: same (weight 13, level 3),
     C23: class number formula of the cubic field x^3-x-1 (L(1,rho) = 2 pi R/sqrt 23, R = log of the plastic number),
     C47a: PARI Hecke L-function of a class-group character of Q(sqrt-47) (bnrinit), identified by lambda(2).
Also: sieved-lattice constants c_P (= ray form truncated) for P = 10, 100, 1000 for all 11 CM/dihedral series,
and direct truncations at p <= 100 for all 15 series (own a_p, tau(p))."""
import os, math, json, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/')
import r1_sums as R

SD = R.SD
PMAX = 30_000_000
import cypari2
P = cypari2.Pari(); P.allocatemem(300 * 10 ** 6)

def primes_to(n):
    fl = R.prime_flags(0, n + 1); p = 1 + np.flatnonzero(fl); return p[p >= 2]

pr = primes_to(PMAX)
fl = R.prime_flags(0, PMAX + 1)
# lambda at split primes for every series of the A, B, C families
lamsq = {}       # key -> (p array of split primes, lambda array)
pA, aA, bA = R.gauss_primes(0, PMAX + 1, fl)
zA = (aA + 1j * bA.astype(float)) / np.sqrt(pA.astype(float))
for m in (4, 8, 12, 16, 20): lamsq['A%d' % m] = (pA, 2 * np.real(zA ** m))
pB, aB, bB = R.eis_primes(0, PMAX + 1, fl)
zB = (aB + 1j * math.sqrt(3) * bB) / np.sqrt(pB.astype(float))
for m in (6, 12): lamsq['B%d' % m] = (pB, 2 * np.real(zB ** m))
odd = pr[pr > 2]
lg23 = R.leg(odd, 23); spl = R.frob_split(odd)
m23 = (lg23 == 1); lamsq['C23'] = (odd[m23], np.where(spl[m23], 2.0, -1.0))
# C31 by the Frobenius test for x^3 + x - 1: substitute (own second routine)
def frob_split_31(p):
    p = p.astype(np.int64)
    def mul(A, B):
        a0, a1, a2 = A; b0, b1, b2 = B
        d0 = (a0 * b0) % p; d1 = ((a0 * b1) % p + (a1 * b0) % p) % p
        d2 = ((a0 * b2) % p + (a1 * b1) % p + (a2 * b0) % p) % p
        d3 = ((a1 * b2) % p + (a2 * b1) % p) % p; d4 = (a2 * b2) % p
        # x^3 = 1 - x, x^4 = x - x^2
        return ((d0 + d3) % p, (d1 - d3 + d4) % p, (d2 - d4) % p)
    def mulx(A):
        a0, a1, a2 = A; return (a2 % p, (a0 - a2) % p, a1 % p)
    Rr = (np.ones_like(p), np.zeros_like(p), np.zeros_like(p))
    for k in range(int(p.max()).bit_length() - 1, -1, -1):
        Rr = mul(Rr, Rr); bit = ((p >> k) & 1).astype(bool); Rx = mulx(Rr); Rr = tuple(np.where(bit, x, y) for x, y in zip(Rx, Rr))
    return (Rr[0] % p == 0) & (Rr[1] % p == 1) & (Rr[2] % p == 0)
lg31 = R.leg(odd, 31); m31 = lg31 == 1; s31 = frob_split_31(odd[m31]); lamsq['C31'] = (odd[m31], np.where(s31, 2.0, -1.0))
lam47 = R.c47_lam(odd); m47 = R.leg(odd, 47) == 1
lamsq['C47a'] = (odd[m47], lam47[m47])
C5a, C5b = R.C5a, R.C5b
lam47b = np.where(np.isclose(lam47, C5a), C5b, np.where(np.isclose(lam47, C5b), C5a, lam47))
lamsq['C47b'] = (odd[m47], lam47b[m47])
LAM2 = {'C23': -1.0, 'C31': -1.0, 'C47a': R.C5a, 'C47b': R.C5b}     # p = 2 splits for D = -23, -31, -47
for k_, l_ in LAM2.items():
    ps_, la_ = lamsq[k_]; lamsq[k_] = (np.concatenate([[2], ps_]), np.concatenate([[l_], la_]))

L1chi = {'A': math.pi / 4, 'B': math.pi / (3 * math.sqrt(3)), 'C23': 3 * math.pi / math.sqrt(23), 'C31': 3 * math.pi / math.sqrt(31), 'C47': 5 * math.pi / math.sqrt(47)}
def l1chi(key): return L1chi[key[0]] if key[0] in 'AB' else L1chi[key[:3]]
D = {'A': 4, 'B': 3}
def disc(key): return D[key[0]] if key[0] in 'AB' else int(key[1:3])

def psi2_lam(key, lam):
    return lam * lam - 2      # lambda of psi^2 at a split prime

out = {'ray': {}, 'direct': {}, 'sieved': {}, 'accel': {}}
checkpoints = [10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, PMAX]
for key, (ps, lam) in lamsq.items():
    o = np.argsort(ps); ps = ps[o]; lam = lam[o]; l2 = psi2_lam(key, lam)
    lr = np.cumsum(np.log1p(-l2 / (ps - 1.0)))
    ray = {}
    for X in checkpoints:
        i = np.searchsorted(ps, X, side='right'); ray[X] = math.exp(lr[i - 1] if i > 0 else 0.0) / l1chi(key)
    out['ray'][key] = ray
    # direct product over all primes (inert 1 + 1/p, ramified 1, split 1 - (lam^2-1)/p); p = 2 for C cases is split (lam(2) given)
    Dd = disc(key)
    allp = pr
    cp = np.ones(len(allp))
    # split mask
    idx = np.searchsorted(allp, ps); cp[idx] = 1 - (lam * lam - 1) / ps
    split = np.zeros(len(allp), bool); split[idx] = True
    ram = (Dd % allp == 0) if key[0] in 'AB' else (allp == Dd)
    if key[0] == 'A': ram = allp == 2
    if key[0] == 'B': ram = allp == 3
    inert = ~split & ~ram
    cp[inert] = 1 + 1.0 / allp[inert]
    lc = np.cumsum(np.log(cp)); dr = {}
    for X in checkpoints:
        i = np.searchsorted(allp, X, side='right'); dr[X] = math.exp(lc[i - 1])
    out['direct'][key] = dr
    # accelerated product (E_p over all p <= PMAX)
    E = np.ones(len(allp))
    l2all = np.zeros(len(allp)); l2all[idx] = l2
    E[split] = (1 - l2all[split] / (allp[split] - 1.0)) / (1 - l2all[split] / allp[split] + 1.0 / allp[split] ** 2)
    E[inert] = 1 / (1 - 1.0 / allp[inert] ** 2)
    E[ram] = 1 / (1 - 1.0 / allp[ram])
    out['accel'][key] = {'prodE': float(np.exp(np.sum(np.log(E)))), 'prodE_1e6': float(np.exp(np.sum(np.log(E[allp <= 10 ** 6]))))}
    print(key, 'ray', {k: round(v, 6) for k, v in ray.items()}, flush=True)
    print(key, 'direct', {k: round(v, 6) for k, v in dr.items()}, flush=True)

# ---- L(1, psi^2) ----------------------------------------------------------------
Lpsi2 = {}
def cm_coeffs_gauss(m2, M):
    """a(n) = (1/4) sum_{a^2+b^2 = n} (a+ib)^m2, n <= M, exact integers (Gaussian integer arithmetic)."""
    a = [0] * (M + 1); r = math.isqrt(M)
    for x in range(-r, r + 1):
        for y in range(-r, r + 1):
            n = x * x + y * y
            if 0 < n <= M:
                z = complex(1, 0); re, im = 1, 0
                for _ in range(m2): re, im = re * x - im * y, re * y + im * x
                a[n] += re
    assert all(v % 4 == 0 for v in a[1:])
    return [v // 4 for v in a[1:]]
def cm_coeffs_eis(m2, M):
    """a(n) = (1/6) sum_{z in Z[w], N z = n} z^m2 (real part; exact), z = x + y w, w^2 = -1 - w, N = x^2 - xy + y^2."""
    a = [0] * (M + 1); r = 2 * math.isqrt(M) + 2
    for x in range(-r, r + 1):
        for y in range(-r, r + 1):
            n = x * x - x * y + y * y
            if 0 < n <= M:
                u, v = 1, 0          # u + v w
                for _ in range(m2):  # (u + v w)(x + y w) = ux + (uy + vx) w + vy w^2 = (ux - vy) + (uy + vx - vy) w
                    u, v = u * x - v * y, u * y + v * x - v * y
                # real part of u + v w = u - v/2 ; imaginary part v sqrt3/2 ; sum over associates is real
                a[n] += 2 * u - v    # store 2 * real part
    assert all(v % 12 == 0 for v in a[1:])
    return [v // 12 for v in a[1:]]  # (1/6) * (2 Re)/2
for key, m, fun, k, N in (('A4', 8, cm_coeffs_gauss, 9, 4), ('B6', 12, cm_coeffs_eis, 13, 3)):
    M = 4000
    co = fun(m, M)
    P('co = %s' % str(co))
    for eps in (1, -1):
        P('Lf = lfuncreate([co, 0, [0,1], %d, %d, %d])' % (k, N, eps))
        fe = int(P('lfuncheckfeq(Lf)'))
        if fe < -30: break
    val = float(P('real(lfun(Lf, %d))' % ((k - 1) // 2 + 1)))
    Lpsi2[key] = {'L1': val, 'feq': fe, 'eps': eps, 'a(1..8)': co[:8]}
    print(key, 'L(1,psi^2) =', val, 'feq', fe, 'eps', eps, co[:8], flush=True)
# C23: class number formula
th = float(P('polrootsreal(x^3-x-1)[1]')); Rg = math.log(th)
Lpsi2['C23'] = {'L1': 2 * math.pi * Rg / math.sqrt(23), 'R': Rg, 'unit': th}
th31 = float(P('polrootsreal(x^3+x-1)[1]')); R31 = abs(math.log(th31))
Lpsi2['C31'] = {'L1': 2 * math.pi * R31 / math.sqrt(31), 'R': R31, 'unit': th31,
                'check_bnf': [float(P('bnfinit(x^3-x-1).reg')), float(P('bnfinit(x^3+x-1).reg')), int(P('bnfinit(x^3-x-1).no')), int(P('bnfinit(x^3+x-1).no'))]}
# C47: Hecke L of class-group characters
P('K = bnfinit(y^2+47); bnr = bnrinit(K, 1); P2 = idealprimedec(K, 2)[1]; d2 = bnrisprincipal(bnr, P2, 0)[1]')
cyc = int(P('bnr.cyc[1]')); d2 = int(P('d2'))
for kk in range(1, cyc):
    l2v = 2 * math.cos(2 * math.pi * kk * d2 / cyc)
    if abs(l2v - C5b) < 1e-9:
        P('Lb = lfuncreate([bnr, [%d]])' % kk); fe = int(P('lfuncheckfeq(Lb)'))
        Lpsi2['C47a'] = {'L1': float(P('real(lfun(Lb, 1))')), 'feq': fe, 'char': kk, 'lam2': l2v}
    if abs(l2v - C5a) < 1e-9:
        P('La = lfuncreate([bnr, [%d]])' % kk); fe = int(P('lfuncheckfeq(La)'))
        Lpsi2['C47b'] = {'L1': float(P('real(lfun(La, 1))')), 'feq': fe, 'char': kk, 'lam2': l2v}
print('Lpsi2', Lpsi2, flush=True)
for key in Lpsi2:
    out['accel'][key]['L1chi'] = l1chi(key); out['accel'][key]['L1psi2'] = Lpsi2[key]['L1']
    out['accel'][key]['c'] = out['accel'][key]['prodE'] / (l1chi(key) * Lpsi2[key]['L1'])
    out['accel'][key]['c_1e6'] = out['accel'][key]['prodE_1e6'] / (l1chi(key) * Lpsi2[key]['L1'])
out['Lpsi2'] = Lpsi2
# ---- sieved-lattice constants (ray form truncated) -----------------------------------
for key in out['ray']:
    out['sieved'][key] = {'lattice(P=0)': 1 / l1chi(key), **{'P=%d' % X: out['ray'][key][X] for X in (10, 100, 1000)}}

# ---- direct truncation at p <= 100 for the non-CM series -----------------------------
def ap_curve(ai, p):
    a1, a2, a3, a4, a6 = ai; cnt = 1
    for x in range(p):
        for y in range(p):
            if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % p == 0: cnt += 1
    return p + 1 - cnt
def tau_list(M):
    c = [0] * (M + 1); c[0] = 1
    for n in range(1, M + 1):          # prod (1 - q^n)^24
        for _ in range(24):
            for j in range(M, n - 1, -1): c[j] -= c[j - n]
    return [0] + c[:M]                 # tau(n) = coefficient of q^n in q prod = c[n-1]
sp100 = [int(x) for x in pr[pr <= 1000]]
trunc = {}
for key, ai, N in (('D11a', (0, -1, 1, -10, -20), 11), ('D37a', (0, 0, 1, -1, 0), 37)):
    cP = {}; prod = 1.0; prodS = 1.0
    for p in sp100:
        if p > 1000: break
        ap = ap_curve(ai, p) if p <= 100 else None
        if p > 100: break
        if p == N: prod *= 1 + 1.0 / p; prodS *= (1 - 1.0 / p ** 3) / (1 - 1.0 / p)
        else:
            lam = ap / math.sqrt(p); prod *= 1 - (lam * lam - 1) / p
            if key == 'D11a':
                s = lam * lam - 2; l1 = s + 1                # Sym^2: Satake {x, 1, 1/x}, s = x + 1/x
                h = [1.0, l1]
                for k in range(2, 4): h.append(l1 * h[k - 1] - l1 * h[k - 2] + (h[k - 3] if k >= 3 else 0.0))
                X = 1.0 / p; prodS *= 1 - (l1 * l1 - 1) * X + h[3] * X * X
    trunc[key] = prod
    if key == 'D11a': trunc['E11a'] = prodS
tl = tau_list(100); prod = 1.0
for p in sp100:
    if p > 100: break
    lam = tl[p] / p ** 5.5; prod *= 1 - (lam * lam - 1) / p
trunc['DD'] = prod
out['direct100_nonCM'] = trunc
print('direct truncations p<=100 (non-CM):', trunc)
json.dump(out, open(SD + '/r3_out.json', 'w'), indent=1, default=lambda o: float(o))
