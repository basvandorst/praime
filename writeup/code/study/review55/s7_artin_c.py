"""c_rho for the Artin factor of a cubic field by the direct product prod_p G_p(1/p) (G_p from L_p(rho) = zeta_{K,p}/zeta_p),
to P = 1e7, with the tail nu/(sqrt P ln P), nu = 2 (Ad rho = chi_d0 + rho), and a log-Cesaro mean.  T1: check c = q/phi(q)."""
import json, math, os, sys
import cypari2
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'; D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json'))
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
# T1 check
bad = []
for k, e in P['L'].items():
    if e['cls'] == 'T1':
        q = e['spec'][1]; ph = int(pari.eulerphi(q))
        if abs(e['c'] - q / ph) > 1e-12: bad.append(k)
print('T1: c = q/phi(q) for all 40:', not bad, bad)
pari("logG(nf, p) = my(fs = [pr.f | pr <- idealprimedec(nf, p)], a, c1, c2, F); a = prod(i=1,#fs, 1 - 'X^fs[i]) / (1 - 'X); c1 = -polcoef(a, 1, 'X); c2 = polcoef(a, 2, 'X); F = 1 - c1^2/p + (c1^2 - c2)*c2/p^2; log(abs(F/(1 - 1/p)))")
pari("runc(pol, PM) = my(nf = nfinit(pol), S = 0., out = List(), chk = [10^4, 10^5, 10^6, 10^7], j = 1, ces = 0., wsum = 0.); forprime(p = 2, PM, S += logG(nf, p); if(p >= PM/10, ces += (S + 2/(sqrt(p)*log(p)))/p; wsum += 1./p); while(j <= #chk && p >= chk[j], listput(out, [chk[j], exp(S), exp(S + 2/(sqrt(chk[j])*log(chk[j])))]); j++)); [Vec(out), exp(ces/wsum)]")
keys = sys.argv[2:] or ['T5_d44', 'T5_d172', 'T5_d331', 'T5_d440', 'T5_d499']
PM = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 6
for k in keys:
    e = P['products'][k]; r = e['rho']
    res = pari('runc(%s, %d)' % (e['spec'][1], PM))
    pts = '  '.join('P=%.0e raw %.5f corr %.5f' % (float(x[0]), float(x[1]), float(x[2])) for x in res[0])
    ces = float(res[1])
    print('%-9s d=%4d reg c_rho = %.5f (c_direct %s) | %s | ces %.5f  ratio ces/reg %.4f' % (k, e['d'], r['c'], r.get('c_direct'), pts, ces, ces / r['c']), flush=True)
