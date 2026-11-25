# (a) T1 (degree 1, Q = q): exact Gonek-lemma prime side S(T) = -i^a q^{-1/2} sum_{n<=qT/2pi} Lambda(n) chi(n) e(-n/q) (theorem-level),
#     and its prime-power part (n = p^k, k >= 2: the secondary term; coherent through p^k mod q).
# (b) T5 d = -484: L(rho) has Q = 484 = 22^2 (the only registered degree-2 case with Q a perfect square): prime-power secondary term
#     of the Artin-law reading R, exact explicit-formula integrals over n = p^k, k >= 2.
import os, sys, math, cmath, json, hashlib, numpy as np, cypari2, sympy
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
BL = os.path.dirname(os.path.abspath(__file__))+'/../../../../../../../writeup/code/study/lowzero/universal/blind/'
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
PB = json.load(open(BL + 'predictions_blind.json'))
res = {'T1': {}, 'T5_d484_R': {}}
num = den = numpp = 0.0
for key, e in sorted(PB['L'].items()):
    if e['cls'] != 'T1': continue
    q, m = int(e['spec'][1]), int(e['spec'][2]); T = float(e['T2']); a = int(round(e['mu'][0]))
    pari('G = znstar(%d, 1)' % q); pari('chi = znconreychar(G, %d)' % m)
    chi = {}
    for n in range(q):
        if math.gcd(n, q) == 1:
            chi[n] = cmath.exp(2j * math.pi * float(pari('chareval(G, chi, %d)' % n)))
    tau = sum(chi[n] * cmath.exp(2j * math.pi * n / q) for n in chi)
    eps_own = tau / ((1j) ** a * math.sqrt(q)); eps_reg = complex(*e['eps'])
    assert abs(eps_own - eps_reg) < 1e-8, (key, eps_own, eps_reg)
    assert (a == 0) == (abs(chi[q - 1] - 1) < 1e-9)
    x = q * T / (2 * math.pi)
    S = 0j; Spp = 0j
    for p in sympy.primerange(2, int(x) + 1):
        if q % p == 0: continue
        pk = p; k = 1
        while pk <= x:
            v = math.log(p) * chi[pk % q] * cmath.exp(-2j * math.pi * pk / q)
            S += v
            if k >= 2: Spp += v
            pk *= p; k += 1
    pref = -(1j ** a) / math.sqrt(q)
    S *= pref; Spp *= pref
    pred = complex(*e['pred'])
    res['T1'][key] = dict(q=q, m=m, pred_registered=[pred.real, pred.imag], exact_prime_side=[S.real, S.imag],
                          prime_power_part=[Spp.real, Spp.imag], ratio_exact_to_registered=[(S / pred).real, (S / pred).imag],
                          sd_along=e['sd']['along'])
    num += (S * pred.conjugate()).real; den += abs(pred) ** 2; numpp += (Spp * pred.conjugate()).real
    print('%-9s eps check ok  pred %8.2f%+8.2fi  exact %8.2f%+8.2fi  (ratio %.4f%+.4fi)  prime-power part %+.2f%+.2fi' % (key, pred.real, pred.imag, S.real, S.imag, (S / pred).real, (S / pred).imag, Spp.real, Spp.imag))
print('T1 pooled: exact prime side / registered law = %.4f ; prime-power part alone shifts the pooled slope by %+.4f' % (num / den, numpp / den))
res['T1_pooled'] = dict(exact_over_law=num / den, prime_power_shift=numpp / den)
# (b) d = -484
e = PB['products']['T5_d484']; pol = e['spec'][1]; d = -484; r = e['rho']; Q = float(r['Q']); mu = r['mu']; T = 300.0
assert abs(Q - 484) < 1e-9
nf = pari.nfinit(pari(pol)); pdisc = abs(int(pari.poldisc(pari(pol))))
tt = np.linspace(1e-9, T, 300001); dt = tt[1] - tt[0]
th = 0.5 * tt * math.log(Q / math.pi ** 2) + sum(np.imag(loggamma((0.5 + mm + 1j * tt) / 2)) for mm in mu)
E2 = np.exp(2j * th)
xmax = Q * (T / (2 * math.pi)) ** 2 * 1.3
P2 = 0j; terms = []
for p in sympy.primerange(2, int(xmax ** 0.5) + 2):
    if pdisc % p == 0:
        Pz = np.poly1d([1.0])
        for pr in pari.idealprimedec(nf, p):
            f_ = int(pr[3]); Pz = Pz * np.poly1d([-1.0] + [0.0] * (f_ - 1) + [1.0])
        qq, rr = np.polydiv(Pz.coeffs, np.array([-1.0, 1.0])); P = list(qq[::-1])
    else:
        nr = len(pari.polrootsmod(pari(pol), p)); P = {3: [1.0, -2.0, 1.0], 0: [1.0, 1.0, 1.0], 1: [1.0, 0.0, -1.0]}[nr]
    roots = [1 / z for z in np.roots(P[::-1])] if len(P) > 1 else []
    pk = p * p; k = 2
    while pk <= xmax:
        Lam = sum(z ** k for z in roots).real * math.log(p)
        if abs(Lam) > 1e-12:
            f = E2 * np.exp(-1j * tt * math.log(pk))
            Iv = ((f.sum() - 0.5 * (f[0] + f[-1])) * dt) / (2 * math.pi)
            v = -Lam * pk ** -0.5 * Iv
            P2 += v; terms.append((p, k, v))
        pk *= p; k += 1
Rreg = e['readings']['R_artin']['pred']
print('T5 d=-484 Artin reading R: registered law %.2f ; prime-power secondary term (k>=2, exact integrals) %+.2f%+.2fi ; sd %.1f'
      % (Rreg, P2.real, P2.imag, e['readings']['R_artin']['sd']))
res['T5_d484_R'] = dict(R_registered=Rreg, prime_power_term=[P2.real, P2.imag], R_amended=[Rreg + P2.real, P2.imag], sd=e['readings']['R_artin']['sd'],
                        note='Q(L(rho)) = 484 = 22^2: n = p^2 terms have kernel e(-2p/22) = e(-p/11), coherent through p mod 11')
json.dump(res, open(OUT + 'a3_level.json', 'w'), indent=1)
