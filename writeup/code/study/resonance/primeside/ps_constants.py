"""ps_constants.py: root numbers, L-values and the constants c = prod_p G_p(1/p) = B/L(1, Ad) for every case of ps_lib.SERIES.

c is computed two ways:
  exact:  c = B / L(1, Ad),  B = prod_p E_p absolutely convergent (E_p = c_p * L_p(1, Ad) = 1 + O(p^-2)), p <= PB,
          L(1, Ad) from PARI:  CM/dihedral: L(1, chi_D) * L(1, psi^2) (psi^2 the CM form / Hecke character of twice the type);
          elliptic curves: L(1, Sym^2 E) (lfunsympow).  Local factors:
             good p:      c_p = 1 - (lam^2 - 1)/p,   E_p = (1 - (lam^2-1)X) / ((1 - X)(1 - (lam^2-2)X + X^2)),  X = 1/p
             ramified p (CM, dihedral; one Satake parameter of modulus 1): c_p = 1, E_p = p/(p-1)
             multiplicative p (11a, 37a): c_p = 1 + 1/p, E_p = c_p * (PARI's local factor of Sym^2 at s = 2)
  direct: the conditionally convergent product prod_{p <= P} c_p (the Hardy-Littlewood singular series), and for the
          CM/dihedral cases its 'ray form' (1/L(1, chi_D)) * prod_{split p <= P} (1 - (lam^2 - 2)/(p - 1)).
Root numbers: PARI lfunrootres; lfuncheckfeq printed.  Output: constants.json, small_lams.npz (lambda(p), p <= 1e5).
usage: ps_constants.py [PB]    (PB default 1e7)"""
import sys, os, json, math, time
import numpy as np
import cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ps_lib as L

PB = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 7
pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6); pari.set_real_precision(38)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'lowzero', 'universal'))
t0 = time.time()

pari('cmfind(N,k,chi,p0,ap0,p1)=my(mf=mfinit([N,k,chi],0),B=mfeigenbasis(mf),r=[]);for(i=1,#B,my(c=mfcoefs(B[i],max(p0,p1)));if(c[p0+1]==ap0&&c[p1+1]==0,r=concat(r,[[mf,B[i]]])));r')

def cm_lfun(N, k, chi, p0, gen, p1):
    """the CM newform of level N, weight k, character chi with a(p0) = 2 Re(gen^(k-1)) and a(p1) = 0; returns Ldata, root number, checkfeq."""
    ap0 = pari('round(2*real((%s)^%d))' % (gen, k - 1))
    r = pari('cmfind(%d,%d,%s,%d,%s,%d)' % (N, k, chi, p0, ap0, p1))
    assert len(r) == 1, (N, k, len(r))
    Lf = pari.lfunmf(r[0][0], r[0][1])
    return Lf, pari.lfunrootres(Lf)[2], int(pari.lfuncheckfeq(Lf)), [int(x) for x in pari.mfcoefs(r[0][1], 12)]

const = {}
# ------------------------------------------------------------------ CM forms
CM = {'A': (4, 'Mod(3,4)', 5, '(2+I)', 3, 2, 'Pi/4'), 'B': (3, 'Mod(2,3)', 7, '(2+sqrt(-3))', 5, 3, 'Pi/(3*sqrt(3))')}
for key, (grp, d, Q, mu, k) in L.SERIES.items():
    if grp not in CM: continue
    N, chi, p0, gen, p1, pram, L1chi = CM[grp]; m = int(key[1:])
    Lf, eps, feq, coefs = cm_lfun(N, k, chi, p0, gen, p1)
    Lg, epsg, feqg, _ = cm_lfun(N, 2 * m + 1, chi, p0, gen, p1)
    L1psi2 = float(pari.real(pari.lfun(Lg, m + 1)))           # analytic s = 1  <->  arithmetic s = (2m+1-1)/2 + 1 = m + 1
    const[key] = dict(group=grp, d=d, Q=Q, mu=mu, k=k, eps=[float(pari.real(eps)), float(pari.imag(eps))], feq=feq, feq_psi2=feqg,
                      coefs_1_12=coefs, L1chi=float(pari(L1chi)), L1psi2=L1psi2, ram=pram,
                      Ad='L(1,chi_D) L(1,psi^2), psi^2 = CM form of weight %d level %d' % (2 * m + 1, N))
    print(key, 'eps', eps, 'feq', feq, feqg, 'L(1,psi^2) =', L1psi2, 'a(1..12)', coefs, '%.0fs' % (time.time() - t0), flush=True)

# ------------------------------------------------------------------ dihedral weight-1 forms: class group characters
for key, pol, chi, D in (('C23', 'y^2+y+6', '[1]', 23), ('C31', 'y^2+y+8', '[1]', 31), ('C47a', 'y^2+y+12', '[1]', 47), ('C47b', 'y^2+y+12', '[2]', 47)):
    bnf = pari.bnfinit(pari(pol), 1); bnr = pari.bnrinit(bnf, 1)
    Lf = pari.lfuncreate(pari('[%s, %s]' % ('bnrinit(bnfinit(%s,1),1)' % pol, chi)))
    an = [complex(pari.real(x), pari.imag(x)) for x in pari.lfunan(Lf, 12)]
    eps = pari.lfunrootres(Lf)[2]; feq = int(pari.lfuncheckfeq(Lf))
    grp, d, Q, mu, k = L.SERIES[key]
    gen = str(pari('bnfinit(%s,1).gen' % pol)); cyc = str(pari('bnfinit(%s,1).cyc' % pol))
    const[key] = dict(group=grp, d=d, Q=Q, mu=mu, k=k, eps=[float(pari.real(eps)), float(pari.imag(eps))], feq=feq,
                      an_1_12=[round(z.real, 6) for z in an], class_gen=gen, cyc=cyc, L1chi=float(pari('lfun(lfuncreate(%d),1)' % -D)), ram=D)
    if D in (23, 31):
        const[key]['L1psi2'] = float(pari.real(pari.lfun(Lf, 1))); const[key]['Ad'] = 'L(1,chi_D) L(1,rho) (psi^2 = psi^-1, Ind psi^-1 = rho)'
    print(key, 'eps', eps, 'feq', feq, 'a(1..12)', [round(z.real, 4) for z in an], 'gen', gen, cyc, flush=True)
# for D = -47: Ad(rho_a) = chi_D + rho_b, Ad(rho_b) = chi_D + rho_a
La = float(pari.real(pari.lfun(pari.lfuncreate(pari('[bnrinit(bnfinit(y^2+y+12,1),1), [1]]')), 1)))
Lb = float(pari.real(pari.lfun(pari.lfuncreate(pari('[bnrinit(bnfinit(y^2+y+12,1),1), [2]]')), 1)))
const['C47a']['L1psi2'] = Lb; const['C47a']['Ad'] = 'L(1,chi_-47) L(1,rho_b)'
const['C47b']['L1psi2'] = La; const['C47b']['Ad'] = 'L(1,chi_-47) L(1,rho_a)'
# which PARI character is our C47a (lam = 2cos(2pi/5) on the class of 2x^2+xy+6y^2, i.e. at p = 2)?
for key in ('C47a', 'C47b'):
    want = L.special_primes(key)[2]
    assert abs(const[key]['an_1_12'][1] - want) < 1e-6, (key, const[key]['an_1_12'][1], want)

# ------------------------------------------------------------------ elliptic curves: c = B / L(1, Sym^2 E)
for key, ai in (('D11a', '[0,-1,1,-10,-20]'), ('D37a', '[0,0,1,-1,0]')):
    E = pari('ellinit(%s)' % ai); N = int(pari.ellglobalred(E)[0]); eps = int(pari.ellrootno(E))
    Ls = pari.lfunsympow(E, 2); L1 = float(pari.real(pari.lfun(Ls, 2)))
    eul = pari('my(E=ellinit(%s)); subst(lfuneuler(lfunsympow(E,2), %d), x, %d^-2)' % (ai, N, N)) if False else None
    loc = float(pari('my(E=ellinit(%s), f=lfuneuler(lfunsympow(E,2), %d)); subst(f, variable(f), 1/%d^2)' % (ai, N, N)))
    grp, d, Q, mu, k = L.SERIES[key]
    const[key] = dict(group=grp, d=d, Q=Q, mu=mu, k=k, eps=[float(eps), 0.0], feq=int(pari.lfuncheckfeq(pari.lfuncreate(E))), L1Ad=L1, bad=N,
                      bad_local_Sym2_at_1=loc, Ad='L(1, Sym^2 E) (lfunsympow, arithmetic s = 2)')
    print(key, 'N', N, 'eps', eps, 'L(1,Sym2) =', L1, 'local factor of Sym2 at the bad prime (s=1):', loc, flush=True)

# ------------------------------------------------------------------ reused constants (round 151 constants files)
for key, f in (('E11a', 'constants_S11a.json'), ('DD', 'constants_D12.json')):
    j = json.load(open(os.path.join(REPO, f))); grp, d, Q, mu, k = L.SERIES[key]
    assert abs(j['Q'] - Q) < 1e-9 and [float(x) for x in j['mu']] == mu
    const[key] = dict(group=grp, d=d, Q=Q, mu=mu, k=k, eps=j['eps'], c_exact=j['c'], c_method='reused: ' + f + ' (' + j['c_method'] + ')')
print('PARI part done %.0fs' % (time.time() - t0), flush=True)

# ------------------------------------------------------------------ lambda(p) for p <= PB from the enumerations; B and direct products
SEG = 2 * 10 ** 6
small = {}
logB = {k: 0.0 for k in const}; logc = {k: 0.0 for k in const}; logray = {k: 0.0 for k in const}
chk = {k: {} for k in const}; CHK = [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
def accumulate(key, p, lam, chi):
    """add primes p (sorted) with eigenvalues lam to the products; ramified/bad handled by caller (chi = 0 there)."""
    X = 1.0 / p.astype(np.float64); a2 = lam * lam
    good = chi != 0
    cp = np.where(good, 1 - (a2 - 1) * X, 1.0); Ep = np.where(good, cp / ((1 - X) * (1 - (a2 - 2) * X + X * X)), 1.0)
    return cp, Ep
allp_cache = {}
def run_products(key, plist, lams):
    o = np.argsort(plist); p = plist[o]; lam = lams[o]
    chi = L.chi_neben(key, p)
    cp, Ep = accumulate(key, p, lam, chi)
    grp = const[key]['group']
    if key[0] in 'ABC':
        pr = const[key]['ram']; ir = np.nonzero(p == pr)[0]
        cp[ir] = 1.0; Ep[ir] = pr / (pr - 1.0)
        ray = np.where((chi == 1), 1 - (lam * lam - 2) / (p - 1.0), 1.0)      # split primes; chi_D part separately
        ray[ir] = 1.0
    else:
        N = const[key]['bad']; ir = np.nonzero(p == N)[0]
        cp[ir] = 1 + 1.0 / N; Ep[ir] = (1 + 1.0 / N) * const[key]['bad_local_Sym2_at_1']
        ray = np.ones_like(cp)
    lc = np.cumsum(np.log(cp)); lr = np.cumsum(np.log(ray))
    for x in CHK:
        i = np.searchsorted(p, x, 'right') - 1
        if i >= 0 and x <= p[-1] + 1: chk[key][x] = [math.exp(lc[i]), math.exp(lr[i])]
    return float(np.sum(np.log(Ep)))

# CM and dihedral cases: enumerate (0, PB]
plists = {k: [] for k in const}; lamlists = {k: [] for k in const}
for Lo in range(0, PB, SEG):
    U = min(PB, Lo + SEG); U += U % 2
    ok = L.sieve_odd(Lo, U)
    for grp, keys in (('A', ['A4', 'A8', 'A12', 'A16', 'A20']), ('B', ['B6', 'B12']), ('C23', ['C23']), ('C31', ['C31']), ('C47', ['C47a', 'C47b'])):
        p, lams, _ = L.lam_series(grp, Lo, U, ok, keys)
        for k in keys: plists[k].append(p); lamlists[k].append(lams[k])
    # all odd primes (for the inert ones in A, B, which carry lam = 0)
    allp = Lo + 1 + 2 * np.nonzero(ok)[0]
    plists['_all'] = plists.get('_all', []) + [allp]
allp = np.concatenate(plists.pop('_all')); lamlists.pop('_all', None)
for key in list(const):
    if key[0] not in 'ABC': continue
    p = np.concatenate(plists[key]); lam = np.concatenate(lamlists[key])
    # add the primes with lam = 0 (inert) and the special primes
    sp = L.special_primes(key)
    if key[0] in 'AB':
        # p, lam are the split primes; inert odd primes: all odd primes not split and not ramified
        q = int(const[key]['Q'] if key[0] == 'B' else 4)
        inert = allp[(allp % q == (3 if q == 4 else 2))]
        p = np.concatenate([p, inert, np.array(list(sp.keys()), np.int64)]); lam = np.concatenate([lam, np.zeros(len(inert)), np.array(list(sp.values()))])
    else:
        sp = {2: sp[2]}          # the odd ramified prime is already in the enumeration (lam set to 1)
        p = np.concatenate([p, np.array([2], np.int64)]); lam = np.concatenate([lam, np.array([sp[2]])])
    logB[key] = run_products(key, p, lam)
    o = np.argsort(p); small[key] = (p[o][p[o] <= 10 ** 5], lam[o][p[o] <= 10 ** 5])
    print(key, 'primes to %.0e: %d (split %d)' % (PB, len(p), int((lam != 0).sum())), '%.0fs' % (time.time() - t0), flush=True)

# elliptic curves: a_p from PARI to PE
PE = min(PB, 10 ** 6)
for key, ai in (('D11a', '[0,-1,1,-10,-20]'), ('D37a', '[0,0,1,-1,0]')):
    v = pari('my(E=ellinit(%s)); [primes([2,%d]), apply(p->ellap(E,p), primes([2,%d]))]' % (ai, PE, PE))
    p = np.array([int(x) for x in v[0]], np.int64); ap = np.array([int(x) for x in v[1]], np.float64)
    lam = ap / np.sqrt(p)
    logB[key] = run_products(key, p, lam)
    small[key] = (p[p <= 10 ** 5], lam[p <= 10 ** 5])
    if key == 'D11a':
        lS = lam * lam - 1; lS[p == 11] = 1.0 / 11; small['E11a'] = (p[p <= 10 ** 5], lS[p <= 10 ** 5])
v = pari('my(c=mfcoefs(mfDelta(), 10^5), P=primes([2,10^5])); [P, apply(p->c[p+1], P)]')
p = np.array([int(x) for x in v[0]], np.int64); tau = np.array([float(int(x)) for x in v[1]])
small['DD'] = (p, tau / p.astype(float) ** 5.5)

# ------------------------------------------------------------------ assemble
for key, cst in const.items():
    if 'c_exact' in cst: continue
    B = math.exp(logB[key])
    L1Ad = cst['L1Ad'] if 'L1Ad' in cst else cst['L1chi'] * cst['L1psi2']
    cst['B'] = B; cst['L1Ad'] = L1Ad; cst['c_exact'] = B / L1Ad
    cst['c_method'] = 'B/L(1,Ad), B over p <= %.0e (%s)' % (PB if key[0] in 'ABC' else PE, cst['Ad'])
    cst['c_direct'] = {('%.0e' % x): v[0] for x, v in chk[key].items()}
    if key[0] in 'ABC': cst['c_rayform'] = {('%.0e' % x): v[1] / cst['L1chi'] for x, v in chk[key].items()}
    print('%-5s eps %+.0f  c = %.8f  (B = %.8f, L(1,Ad) = %.8f)   direct products: %s   ray form: %s' % (
        key, cst['eps'][0], cst['c_exact'], B, L1Ad, ' '.join('%s:%.5f' % kv for kv in cst['c_direct'].items()),
        ' '.join('%s:%.5f' % kv for kv in cst.get('c_rayform', {}).items())), flush=True)
json.dump(const, open(os.path.join(HERE, 'constants.json'), 'w'), indent=1)
np.savez(os.path.join(HERE, 'small_lams.npz'), **{k + '_p': v[0] for k, v in small.items()}, **{k + '_lam': v[1] for k, v in small.items()})
print('done %.0fs' % (time.time() - t0))
