"""blind test of the one-point law: every registered reading with its decision.
usage: bt_test.py            (predictions_blind.json, zeros_blind/)   -- refuses a prediction file whose sha256 != FROZEN
       bt_test.py --dry      (dry/predictions_dry.json, dry/zeros/)   -- refuses unless sha256 == DRY
Excluded: keys without a COMPLETE zero file (listed). Pairs read the round-144 zero lists (repository, read-only).
Writes results_blind.json (results_dry.json with --dry)."""
import sys, os, json, math, cmath, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from ulib import theta
FROZEN = 'c8a77047709d6e8c4fe5dbfb1fa595961dc3ef25ef9d2082ea7245eece3bf50f'
DRY = 'f48453daaaf9725f76b994345a48c12c0e6f69d4f09eca88dc16d45537930456'
R144 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..', 'data', 'study', 'r144', 'fresh_zeros.npz')
dry = '--dry' in sys.argv
pf, zdir, want, outf = ('dry/predictions_dry.json', 'dry/zeros', DRY, 'results_dry.json') if dry else ('predictions_blind.json', 'zeros_blind', FROZEN, 'results_blind.json')
raw = open(pf, 'rb').read(); h = hashlib.sha256(raw).hexdigest()
if h != want:
    sys.exit('REFUSED: sha256(%s) = %s, registered %s' % (pf, h, want))
P = json.loads(raw); OUT = dict(predictions=pf, sha256=h, readings={}, excluded=[])
print('prediction file %s sha256 %s (matches the registered hash)' % (pf, h))
def load(key):
    fn = os.path.join(zdir, key + '.npz')
    if not os.path.exists(fn): return None, dict(status='MISSING')
    z = np.load(fn); info = json.loads(str(z['info']))
    return (z['g'] if info.get('status') == 'COMPLETE' else None), info
def dec(ok): return 'PASS' if ok else 'FAIL'
LINES = {'2': 2.0, '1/2': 0.5, '3': 3.0, '1/3': 1 / 3}
# ---------------------------------------------------------------- primitive L
per = {}
for key, e in P['L'].items():
    g, info = load(key)
    if g is None: OUT['excluded'].append([key, info.get('status'), info.get('reason', '')]); continue
    T2 = e['T2']; g = g[(g > 0) & (g <= T2)]
    ph = np.exp(2j * theta(g, e['Q'], e['mu'])); w = np.sin(math.pi * g / T2) ** 2
    G = complex(ph.sum()); Gh = complex((w * ph).sum())
    Gl = {nm: complex((ph * np.exp(-1j * g * math.log(r))).sum()) for nm, r in LINES.items()}
    pred = complex(*e['pred'])
    per[key] = dict(cls=e['cls'], G=G, Gh=Gh, Gl=Gl, pred=pred, n=len(g), info=info)
    print('%-30s %-4s n=%4d  pred %+8.1f%+7.1fi  G %+8.1f%+7.1fi  slope %+.3f ± %.3f  Hann %+.3f  [%s, drift %+.3f]' % (
        key, e['cls'], len(g), pred.real, pred.imag, G.real, G.imag, (G / pred).real if abs(pred) > 0 else float('nan'),
        e['sd']['along'] / abs(pred) if abs(pred) > 0 else float('nan'), (Gh / (pred / 2)).real if abs(pred) > 0 else float('nan'), info['status'], info['drift']))
def pooled(keys, which='r1'):
    num = var = 0.0; numc = 0j; vim = 0.0
    for k in keys:
        e = P['L'][k]; r = per[k]
        if which == 'r1':
            pr = complex(*e['pred']); G = r['G']; sa, sc = e['sd']['along'], e['sd']['across']
            num += abs(pr) ** 2; numc += G * pr.conjugate(); var += abs(pr) ** 2 * sa ** 2; vim += abs(pr) ** 2 * sc ** 2
        elif which == 'hann':
            pr = complex(*e['pred_hann']); G = r['Gh']; sa = e['sd']['along_hann']
            num += abs(pr) ** 2; numc += G * pr.conjugate(); var += abs(pr) ** 2 * sa ** 2
        else:
            for nm in LINES:
                pr = complex(*e['lines'][nm]['pred']); G = r['Gl'][nm]; sa = e['sd']['line_along']
                num += abs(pr) ** 2; numc += G * pr.conjugate(); var += abs(pr) ** 2 * sa ** 2
    if num == 0: return None
    return dict(n=len(keys), slope=numc.real / num, im=numc.imag / num, se=math.sqrt(var) / num, se_im=math.sqrt(vim) / num if vim else None)
fams = sorted(set(e['cls'] for e in P['L'].values()))
print('\nD1 family slopes (noise-model se; PASS iff |b - 1| <= 3 se):')
for f in fams:
    ks = [k for k in per if per[k]['cls'] == f]
    if not ks: print('  %s: no complete L' % f); continue
    r = pooled(ks); rh = pooled(ks, 'hann')
    OUT['readings']['D1_' + f] = dict(r, hann=rh, decision=dec(abs(r['slope'] - 1) <= 3 * r['se']), falsified=abs(r['slope'] - 1) > max(4 * r['se'], 0.10))
    print('  %-4s n=%2d  b = %.4f ± %.4f  (Hann %.4f ± %.4f)  %s' % (f, r['n'], r['slope'], r['se'], rh['slope'], rh['se'], OUT['readings']['D1_' + f]['decision']))
allk = list(per)
if allk:
    r = pooled(allk)
    OUT['readings']['D2_all'] = dict(r, decision=dec(abs(r['slope'] - 1) <= 3 * r['se']), falsified=abs(r['slope'] - 1) > max(4 * r['se'], 0.05))
    OUT['readings']['D3_H0'] = dict(z=r['slope'] / r['se'], decision='H0 rejected' if r['slope'] / r['se'] >= 5 else 'H0 NOT rejected')
    print('D2 all primitive n=%d: b = %.4f ± %.4f  %s ;  D3: b/se = %.1f -> %s' % (r['n'], r['slope'], r['se'], OUT['readings']['D2_all']['decision'], r['slope'] / r['se'], OUT['readings']['D3_H0']['decision']))
    rl = pooled(allk, 'lines')
    OUT['readings']['D4_lines'] = dict(rl, decision=dec(abs(rl['slope'] - 1) <= 3 * rl['se']))
    print('D4 lines r = 2, 1/2, 3, 1/3: b = %.4f ± %.4f  %s' % (rl['slope'], rl['se'], OUT['readings']['D4_lines']['decision']))
    for nm in LINES:
        num = var = 0.0
        for k in allk:
            e = P['L'][k]; pr = complex(*e['lines'][nm]['pred']); G = per[k]['Gl'][nm]
            num += (G * pr.conjugate()).real; var += abs(pr) ** 2 * e['sd']['line_along'] ** 2
        den = sum(abs(complex(*P['L'][k]['lines'][nm]['pred'])) ** 2 for k in allk)
        print('    line %-3s b = %.3f ± %.3f' % (nm, num / den, math.sqrt(var) / den))
    # D5 amplitude follows c
    X = []; Y = []; Wt = []
    for k in allk:
        e = P['L'][k]; eps = complex(*e['eps']); unit = e['T2'] / (2 * math.pi)
        y = (per[k]['G'] * (-eps).conjugate()).real / abs(eps) / unit; s = e['sd']['along'] / unit
        X.append([1.0, e['c'] - 1.0]); Y.append(y); Wt.append(1 / s ** 2)
    X = np.array(X); Y = np.array(Y); Wt = np.array(Wt)
    A = X.T @ (Wt[:, None] * X); bvec = np.linalg.solve(A, X.T @ (Wt * Y)); cov = np.linalg.inv(A)
    al, be = bvec; sal, sbe = np.sqrt(np.diag(cov))
    OUT['readings']['D5_c'] = dict(alpha=al, se_alpha=sal, beta=be, se_beta=sbe, decision=dec(abs(be - 1) <= 3 * sbe), H_c1='rejected' if be / sbe >= 3 else 'NOT rejected')
    print('D5 amplitude follows c: alpha = %.3f ± %.3f, beta = %.3f ± %.3f  %s ; H-c1 %s' % (al, sal, be, sbe, OUT['readings']['D5_c']['decision'], OUT['readings']['D5_c']['H_c1']))
    neg = [k for k in allk if P['L'][k]['eps'][0] < -0.5 and abs(P['L'][k]['eps'][1]) < 1e-9]
    if neg:
        r = pooled(neg)
        OUT['readings']['D6_eps_minus'] = dict(r, decision=dec(abs(r['slope'] - 1) <= 3 * r['se']), H_noroot_value=-1.0)
        print('D6 eps = -1 subset n=%d: b = %.3f ± %.3f  %s  (H-noroot: -1)' % (r['n'], r['slope'], r['se'], OUT['readings']['D6_eps_minus']['decision']))
    cpx = [k for k in allk if abs(P['L'][k]['eps'][1]) > 1e-9]
    if cpx:
        r = pooled(cpx)
        num = sum(abs(complex(*P['L'][k]['pred'])) ** 2 for k in cpx)
        hc = sum(complex(*P['L'][k]['alt']['H_conj']) * complex(*P['L'][k]['pred']).conjugate() for k in cpx) / num
        hn = sum(complex(*P['L'][k]['alt']['H_noroot']) * complex(*P['L'][k]['pred']).conjugate() for k in cpx) / num
        OUT['readings']['D7_complex_eps'] = dict(r, H_conj_s=[hc.real, hc.imag], H_noroot_s=[hn.real, hn.imag], decision=dec(abs(r['im']) <= 3 * r['se_im']))
        print('D7 complex eps n=%d: s = %.4f %+.4fi (se %.4f, %.4f)  %s ;  H-conj would give %.3f%+.3fi, H-noroot %.3f%+.3fi' % (
            r['n'], r['slope'], r['im'], r['se'], r['se_im'], OUT['readings']['D7_complex_eps']['decision'], hc.real, hc.imag, hn.real, hn.imag))
# ---------------------------------------------------------------- products
gz, zinfo = load('Zeta')
t5 = {}
for key, e in P['products'].items():
    g, info = load(key)
    if g is None: OUT['excluded'].append([key, info.get('status'), info.get('reason', '')]); continue
    T2 = e['T2']; g = g[(g > 0) & (g <= T2)]
    thK = theta(g, e['Q'], e['mu'])
    if e['cls'] == 'T5':
        if gz is None: OUT['excluded'].append([key, 'no Zeta zeros']); continue
        zz = gz[gz <= T2 + 0.5]
        isz = np.array([np.min(np.abs(zz - x)) < 1e-6 for x in g])
        nzeta = int(np.sum(gz <= T2))
        r = e['rho']; gr = g[~isz]; gzz = g[isz]
        w = lambda x: np.sin(math.pi * x / T2) ** 2
        M1 = complex(np.exp(2j * thK).sum()); M1h = complex((w(g) * np.exp(2j * thK)).sum())
        M2 = complex(np.exp(2j * theta(gr, 1.0, [0.0])).sum()); M2h = complex((w(gr) * np.exp(2j * theta(gr, 1.0, [0.0]))).sum())
        M3 = complex(np.exp(2j * theta(gzz, e['Q'], e['mu'])).sum())
        Ra = complex(np.exp(2j * theta(gr, r['Q'], r['mu'])).sum())
        t5[key] = dict(M1=M1, M1h=M1h, M2=M2, M2h=M2h, M3=M3, R=Ra, n_zeta_matched=int(isz.sum()), n_zeta=nzeta)
        rd = e['readings']
        print('%-12s d=%5d zeta matched %d/%d  M1 %+6.1f (mech %+6.1f)  M2 %+6.1f (mech %+6.1f)  M3 %+6.1f  Artin %+7.1f (law %+7.1f)' % (
            key, e['d'], int(isz.sum()), nzeta, M1.real, rd['M1']['pred'], M2.real, rd['M2']['pred'], M3.real, Ra.real, rd['R_artin']['pred']))
    else:
        half = complex(np.exp(1j * thK).sum()); prod = complex(np.exp(2j * thK).sum())
        t5[key] = dict(half=half, prod=prod)
        print('%-12s half phase %+.1f%+.1fi (pred %+.1f)  product phase %+.1f%+.1fi (pred 0, sd %.1f)' % (key, half.real, half.imag, e['readings']['half_phase']['pred'], prod.real, prod.imag, e['readings']['product_phase']['sd']))
def slope_on(keys, rd, meas):
    num = sum(t5[k][meas].real * P['products'][k]['readings'][rd]['pred'] for k in keys)
    den = sum(P['products'][k]['readings'][rd]['pred'] ** 2 for k in keys)
    var = sum(P['products'][k]['readings'][rd]['pred'] ** 2 * P['products'][k]['readings'][rd]['sd'] ** 2 for k in keys)
    return (num / den, math.sqrt(var) / den) if den > 0 else (float('nan'), float('nan'))
k5 = [k for k in t5 if P['products'][k]['cls'] == 'T5']
if k5:
    for rd, meas, D in (('M1', 'M1', 'D8_M1'), ('M2', 'M2', 'D9_M2'), ('R_artin', 'R', 'D10_artin')):
        b, s = slope_on(k5, rd, meas)
        OUT['readings'][D] = dict(n=len(k5), slope=b, se=s, decision=dec(abs(b - 1) <= 3 * s), H_prod0='rejected' if b / s >= 3 else 'NOT rejected')
        print('%s n=%d: slope on the prediction %.3f ± %.3f  %s ; zero-hypothesis %s' % (D, len(k5), b, s, OUT['readings'][D]['decision'], OUT['readings'][D]['H_prod0']))
    for k in k5:                      # D10b: primitive Artin L with c_rho = 0 exactly (F_p(1/p) = 0, p = 3 splits completely)
        e = P['products'][k]
        if abs(e['rho']['c']) < 1e-9:
            sdz = e['readings']['R_artin']['sd']; Rv = t5[k]['R']
            OUT['readings']['D10b_c0_' + k] = dict(R=[Rv.real, Rv.imag], sd=sdz, decision=dec(abs(Rv.real) <= 3 * sdz), H_c1=-e['T2'] / (2 * math.pi))
            print('D10b %s (c_rho = 0): Artin sum %+.1f%+.1fi against 0 ± %.1f  %s  (a c = 1 law would give %+.1f)' % (k, Rv.real, Rv.imag, sdz, OUT['readings']['D10b_c0_' + k]['decision'], -e['T2'] / (2 * math.pi)))
    m3 = np.array([t5[k]['M3'].real for k in k5]); sd3 = np.mean([P['products'][k]['readings']['M3']['sd'] for k in k5]) / math.sqrt(len(k5))
    OUT['readings']['D9_M3'] = dict(mean=float(m3.mean()), se=sd3, decision=dec(abs(m3.mean()) <= 3 * sd3))
    print('D9 M3 mean %.2f ± %.2f  %s' % (m3.mean(), sd3, OUT['readings']['D9_M3']['decision']))
kg = [k for k in t5 if P['products'][k]['cls'] == 'T5g']
if kg:
    num = sum(t5[k]['half'].real * P['products'][k]['readings']['half_phase']['pred'] for k in kg)
    den = sum(P['products'][k]['readings']['half_phase']['pred'] ** 2 for k in kg)
    var = sum(P['products'][k]['readings']['half_phase']['pred'] ** 2 * P['products'][k]['readings']['half_phase']['sd'] ** 2 for k in kg)
    b, s = num / den, math.sqrt(var) / den
    pz = [abs(t5[k]['prod'].real) <= 3 * P['products'][k]['readings']['product_phase']['sd'] for k in kg]
    OUT['readings']['D11_T5g'] = dict(half_slope=b, se=s, decision=dec(abs(b - 1) <= 3 * s and all(pz)))
    print('D11 GL2-type products: half-phase slope %.3f ± %.3f, product phase within 3 sd: %s  %s' % (b, s, pz, OUT['readings']['D11_T5g']['decision']))
# ---------------------------------------------------------------- pairs (round-144 zero lists)
if P['pairs']:
    z = np.load(R144); vals = []
    for pr in P['pairs']:
        T2 = pr['T2']
        g1 = z[pr['keys'][0]]; g2 = z[pr['keys'][1]]; g1 = g1[(g1 > 0) & (g1 <= T2)]; g2 = g2[(g2 > 0) & (g2 <= T2)]
        th = lambda t, N: theta(t, float(N), [0.5, 1.5])
        gu = np.concatenate([g1, g2])
        S = complex(np.exp(2j * (th(gu, pr['N'][0]) + th(gu, pr['N'][1]))).sum())
        vals.append([S, math.sqrt(len(gu))])
        print('pair %s + %s: product-phase sum %+.1f%+.1fi (pred 0, sd %.1f)' % (pr['keys'][0], pr['keys'][1], S.real, S.imag, math.sqrt(len(gu))))
    mre = np.mean([v[0].real for v in vals]); mim = np.mean([v[0].imag for v in vals]); sdm = np.mean([v[1] for v in vals]) / math.sqrt(len(vals))
    OUT['readings']['D12_pairs'] = dict(mean=[mre, mim], se=sdm, decision=dec(abs(mre) <= 3 * sdm and abs(mim) <= 3 * sdm))
    print('D12 pairs: mean %.2f%+.2fi ± %.2f  %s' % (mre, mim, sdm, OUT['readings']['D12_pairs']['decision']))
fal = [k for k, v in OUT['readings'].items() if isinstance(v, dict) and v.get('falsified')]
print('\nexcluded:', OUT['excluded'] if OUT['excluded'] else 'none')
print('FALSIFIER (all-primitive slope outside 1 ± max(4 se, 0.05), or a family outside 1 ± max(4 se_F, 0.10)):', 'TRIGGERED by ' + ', '.join(fal) if fal else 'not triggered')
print('strict statistical decisions (3 se, no allowance) failing:', [k for k, v in OUT['readings'].items() if isinstance(v, dict) and v.get('decision') == 'FAIL'] or 'none')
def cj(x):
    if isinstance(x, complex): return [x.real, x.imag]
    if isinstance(x, dict): return {k: cj(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [cj(v) for v in x]
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    return x
OUT['per_L'] = cj({k: dict(G=v['G'], Gh=v['Gh'], Gl=v['Gl'], n=v['n']) for k, v in per.items()}); OUT['products'] = cj(t5)
OUT['readings'] = cj(OUT['readings'])
json.dump(OUT, open(outf, 'w'), indent=1)
