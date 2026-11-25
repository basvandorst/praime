import os, json, math, numpy as np
from scipy.special import loggamma
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
def theta(t, Q, mu):
    t = np.asarray(t, float); d = len(mu)
    return 0.5 * t * math.log(Q / math.pi ** d) + sum(np.imag(loggamma((0.5 + m + 1j * t) / 2)) for m in mu)
a = json.load(open(OUT + 'r7f_n13_T400.json'))['res']; b = json.load(open(OUT + 'r7g_n13_400_600_fixed.json'))['res']
sets = {}
for k in ('N13', 'N13d'):
    g = np.sort(np.concatenate([a[k]['zeros'], b[k]['zeros']])); e = complex(*a[k]['eps'])
    sets[k] = dict(g=g, T=600.0, Q=13.0, mu=[0.5, 1.5], eps=e, c=1.0215613, A1=1.1210163, P2=0j)
for tag, fn, Q, mu, c, A1, P2 in (('F24a', 'r10_f24a_T400.json', 1.0, [11.5, 12.5], -0.0857724, 1.8846368, -1.23 - 5.90j),
                                  ('F24b', 'r10g_F24b_T400.json', 1.0, [11.5, 12.5], 0.5270993, 1.5755353, -0.59 - 6.21j),
                                  ('W29a', 'r10g_W29a_T400.json', 29.0, [0.5, 1.5], 0.0396719, 2.6872059, 0j),
                                  ('W29b', 'r10g_W29b_T400.json', 29.0, [0.5, 1.5], 2.3228657, 0.6855927, 0j)):
    r = json.load(open(OUT + fn)); sets[tag] = dict(g=np.array(r['zeros']), T=float(r['T']), Q=Q, mu=mu, eps=1 + 0j, c=c, A1=A1, P2=P2)
rows = []
for k, s in sets.items():
    g = s['g']; T = s['T']; th = theta(g, s['Q'], s['mu']); e = np.exp(2j * th)
    n = len(g); M = theta(T, s['Q'], s['mu']) / math.pi
    W = T / (2 * math.pi); law = -s['eps'] * s['c'] * W
    G = e.sum()
    # block errors of the complex sum (K-averaged variance of Re and Im)
    vr = []; vi = []
    for K in (4, 6, 8, 12, 20, 30):
        edges = np.linspace(0, T, K + 1)
        Gk = np.array([e[(g > edges[j]) & (g <= edges[j + 1])].sum() for j in range(K)])
        vr.append(Gk.real.var(ddof=1) * K); vi.append(Gk.imag.var(ddof=1) * K)
    sr, si = math.sqrt(np.mean(vr)), math.sqrt(np.mean(vi))
    alts = {'law': law, 'law+P2': law + s['P2'], 'c=1': -s['eps'] * W, 'no B (1/L(1,Ad))': -s['eps'] / s['A1'] * W,
            'no root number': -s['c'] * W, 'conj eps': -np.conj(s['eps']) * s['c'] * W, 'H0': 0j}
    print('%-5s T=%.0f zeros %d (n - theta(T)/pi = %+.2f) G = %7.2f%+7.2fi  se(Re) %.1f se(Im) %.1f' % (k, T, n, n - M, G.real, G.imag, sr, si))
    for name, v in alts.items():
        chi = ((G - v).real / sr) ** 2 + ((G - v).imag / si) ** 2
        print('      %-18s %8.2f%+8.2fi   chi2 = %7.2f' % (name, v.real, v.imag, chi))
    rows.append((k, G, alts, sr, si))
# combined phase test for the nebentypus pair: antisymmetric imaginary part
Gf = [r for r in rows if r[0] == 'N13'][0]; Gd = [r for r in rows if r[0] == 'N13d'][0]
anti = (Gf[1].imag - Gd[1].imag) / 2; se_anti = math.sqrt(Gf[4] ** 2 + Gd[4] ** 2) / 2
sym = (Gf[1].real + Gd[1].real) / 2; se_sym = math.sqrt(Gf[3] ** 2 + Gd[3] ** 2) / 2
la = (Gf[2]['law'].imag - Gd[2]['law'].imag) / 2; ls = (Gf[2]['law'].real + Gd[2]['law'].real) / 2
print('nebentypus pair, T = 600: (Im G_f - Im G_fbar)/2 = %.2f +- %.2f ; law %+.2f, no root number 0, conj(eps) %+.2f' % (anti, se_anti, la, -la))
print('                          (Re G_f + Re G_fbar)/2 = %.2f +- %.2f ; law %.2f (slope %.3f +- %.3f)' % (sym, se_sym, ls, sym / ls, se_sym / abs(ls)))
for k in ('N13', 'N13d'):
    r = [x for x in rows if x[0] == k][0]
    print('   %s: arg G = %.1f deg (law %.1f, no root 180, conj %.1f), angular se ~ %.1f deg' % (k, np.degrees(np.angle(r[1])), np.degrees(np.angle(r[2]['law'])), 180.0, np.degrees(np.angle(r[2]['conj eps'])), np.degrees(r[4] / abs(r[1]))))
