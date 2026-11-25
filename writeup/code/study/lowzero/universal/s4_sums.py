"""the zero sums against the frozen predictions. usage: s4_sums.py <sha256 of predictions.json>
Reads predictions.json (refuses if the hash differs), zeros/zeros_<key>.npz; writes results.json, results_table.md."""
import sys, os, json, math, cmath, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ulib import theta
raw = open('predictions.json').read()
h = hashlib.sha256(raw.encode()).hexdigest()
assert h == sys.argv[1], 'predictions.json changed: %s' % h
P = json.loads(raw)
def load(key):
    fn = 'zeros/zeros_%s.npz' % key
    if not os.path.exists(fn): return None, None
    z = np.load(fn); return z['g'], json.loads(str(z['info']))
NB = 6
res = {}; rows = []
for key, e in P['Lfunctions'].items():
    g, info = load(key)
    if g is None: continue
    T2 = e['T2']; Q = e['Q']; mu = e['mu']
    g = g[(g > 0) & (g <= T2)]
    ph = np.exp(2j * theta(g, Q, mu))
    pred = complex(*e['pred_sharp']); eps = complex(*e['eps']); c = e['c']
    G = ph.sum(); Gh = (np.sin(math.pi * g / T2) ** 2 * ph).sum()
    edges = np.linspace(0, T2, NB + 1); bl = np.array([ph[(g > edges[i]) & (g <= edges[i + 1])].sum() for i in range(NB)])
    unitb = (T2 / NB) / (2 * math.pi)
    r = dict(key=key, family=e['family'], n=len(g), G=[G.real, G.imag], G_hann=[Gh.real, Gh.imag], blocks=[[b.real, b.imag] for b in bl],
             completeness=info)
    if c != 0:
        predb = -eps * c * unitb
        s_k = np.array([(b / predb).real for b in bl])
        b = (G / pred).real; bh = (Gh / (pred / 2)).real
        se_blk = float(s_k.std(ddof=1) / math.sqrt(NB)); se_heur = e['sd_heur_sharp'] / abs(pred); se_heur_h = e['sd_heur_hann'] / abs(pred / 2)
        im = (G / pred).imag
        ang = math.degrees(cmath.phase(G / pred))
        r.update(slope=b, slope_hann=bh, se_block=se_blk, se_heur=se_heur, se_heur_hann=se_heur_h, im_ratio=im, angle_deg=ang,
                 angle_se_deg=math.degrees(e['sd_heur_sharp'] / abs(pred)),
                 alt_noroot_slope=(G / complex(*e['alt_noroot'])).real, alt_conj_slope=(G / complex(*e['alt_conj'])).real if abs(complex(*e['alt_conj'])) > 0 else None,
                 PASS=bool(abs(b - 1) <= 3 * se_blk))
    else:
        unit = T2 / (2 * math.pi)
        chat = (-G / (eps * unit)).real; chat_k = np.array([(-bb / (eps * unitb)).real for bb in bl])
        se_blk = float(chat_k.std(ddof=1) / math.sqrt(NB)); se_heur = e['sd_heur_sharp'] / unit
        r.update(c_hat=chat, c_hat_im=(-G / (eps * unit)).imag, se_block=se_blk, se_heur=se_heur, c_naive=e['c_direct_end_1e5'],
                 PASS=bool(abs(chat) <= 3 * se_blk))
    # POST HOC (not registered): empirical noise of the Hann sum from 40 control frequencies xi in +-[0.06, 1.4]
    xis = np.concatenate([np.linspace(0.06, 1.4, 20), -np.linspace(0.06, 1.4, 20)])
    wh = np.sin(math.pi * g / T2) ** 2
    ctrl = np.array([(wh * ph * np.exp(-1j * g * x)).sum() for x in xis])
    r['ctrl_sd_hann'] = float(np.sqrt(np.mean(ctrl.real ** 2)))
    r['ctrl_sd_hann_im'] = float(np.sqrt(np.mean(ctrl.imag ** 2)))
    if c != 0:
        r['z_hann_ctrl'] = float(((Gh - pred / 2) / r['ctrl_sd_hann']).real)
        r['slope_hann_se_ctrl'] = float(r['ctrl_sd_hann'] / abs(pred / 2))
    else:
        r['c_hat_hann'] = float((-Gh / (eps * T2 / (4 * math.pi))).real)
        r['c_hat_hann_se_ctrl'] = float(r['ctrl_sd_hann'] / (T2 / (4 * math.pi)))
    # lines
    L = {}
    for name, rr in (('2', 2.0), ('1/2', 0.5), ('3', 3.0), ('1/3', 1 / 3)):
        Gr = (ph * np.exp(-1j * g * math.log(rr))).sum()
        L[name] = dict(G=[Gr.real, Gr.imag], pred=e['lines'][name]['sharp'])
    r['lines'] = L
    res[key] = r
# family pooling
fams = {}
for key, r in res.items():
    if 'slope' in r: fams.setdefault(r['family'], []).append(key)
pool = {}
for fam, keys in fams.items():
    num = den = vb = vh = 0.0
    for k in keys:
        pr = complex(*P['Lfunctions'][k]['pred_sharp']); G = complex(*res[k]['G'])
        num += (G * pr.conjugate()).real; den += abs(pr) ** 2
        vb += abs(pr) ** 4 * res[k]['se_block'] ** 2; vh += abs(pr) ** 4 * res[k]['se_heur'] ** 2
    pool[fam] = dict(keys=keys, slope=num / den, se_block=math.sqrt(vb) / den, se_heur=math.sqrt(vh) / den)
num = den = vb = vh = 0.0
for k, r in res.items():
    if 'slope' not in r: continue
    pr = complex(*P['Lfunctions'][k]['pred_sharp']); G = complex(*r['G'])
    num += (G * pr.conjugate()).real; den += abs(pr) ** 2; vb += abs(pr) ** 4 * r['se_block'] ** 2; vh += abs(pr) ** 4 * r['se_heur'] ** 2
pool['ALL primitive'] = dict(slope=num / den, se_block=math.sqrt(vb) / den, se_heur=math.sqrt(vh) / den, n=sum(1 for r in res.values() if 'slope' in r))
# lines pooled (primitive L only)
ln = {}
for name in ('2', '1/2', '3', '1/3', 'all'):
    num = den = var = 0.0
    for k, r in res.items():
        if 'slope' not in r: continue
        sd = P['Lfunctions'][k]['sd_heur_sharp']
        for nm in (('2', '1/2', '3', '1/3') if name == 'all' else (name,)):
            pr = complex(*r['lines'][nm]['pred']); G = complex(*r['lines'][nm]['G'])
            num += (G * pr.conjugate()).real; den += abs(pr) ** 2; var += abs(pr) ** 2 * sd ** 2
    ln[name] = dict(slope=num / den, se_heur=math.sqrt(var) / den)
pool['lines'] = ln
# cross structure
cross = {}
gz, _ = load('Zeta')
for Z, R in (('Z23', 'R23'), ('Z31', 'R31')):
    gK, _ = load(Z); gR, _ = load(R)
    if gK is None or gz is None: continue
    eK = P['Lfunctions'][Z]; T2 = eK['T2']
    gK = gK[gK <= T2]; gzz = gz[gz <= T2 + 0.5]
    isz = np.array([np.min(np.abs(gzz - x)) < 1e-6 for x in gK])
    gr_from_K = gK[~isz]
    rinfo = None
    if gR is not None:
        gRR = gR[gR <= T2]
        matched = sum(np.min(np.abs(gRR - x)) < 1e-6 for x in gr_from_K)
        rinfo = dict(rho_zeros_in_K=len(gr_from_K), rho_zeros_direct=len(gRR), matched=int(matched))
    eR = P['Lfunctions'][R]
    thK = lambda t: theta(t, eK['Q'], eK['mu']); thz = lambda t: theta(t, 1.0, [0.0]); thR = lambda t: theta(t, eR['Q'], eR['mu'])
    S = lambda gg, th: complex(np.exp(2j * th(gg)).sum())
    zz = gK[isz]
    cross[Z] = {'counts': dict(all=len(gK), zeta=int(isz.sum()), rho=int((~isz).sum()), zeta_direct=int(np.sum(gz <= T2)), rho_check=rinfo),
                'all zeros, phase theta_K = theta_zeta + theta_rho': S(gK, thK), 'zeta zeros, phase theta_zeta (degree-1 theorem)': S(zz, thz),
                'rho zeros, phase theta_rho (degree-2 law)': S(gr_from_K, thR), 'zeta zeros, phase theta_K': S(zz, thK),
                'rho zeros, phase theta_K': S(gr_from_K, thK), 'all zeros, phase theta_zeta': S(gK, thz), 'all zeros, phase theta_rho': S(gK, thR)}
gJ, _ = load('J13'); gf, _ = load('N13'); gfb, _ = load('N13d')
if gJ is not None and gf is not None and gfb is not None:
    eJ = P['Lfunctions']['J13']; T2 = eJ['T2']; eN = P['Lfunctions']['N13']
    gJ = gJ[gJ <= T2]; gf = gf[gf <= T2]; gfb = gfb[gfb <= T2]
    un = np.sort(np.concatenate([gf, gfb]))
    match = sum(np.min(np.abs(un - x)) < 1e-6 for x in gJ)
    thJ = lambda t: theta(t, eJ['Q'], eJ['mu']); thf = lambda t: theta(t, eN['Q'], eN['mu'])
    S = lambda gg, th: complex(np.exp(2j * th(gg)).sum())
    cross['J13'] = {'counts': dict(J=len(gJ), f=len(gf), fbar=len(gfb), J_matched_in_union=int(match)),
                    'all zeros, phase theta_J = 2 theta_f': S(gJ, thJ), 'all zeros, half phase theta_f (union of the f and f-bar laws)': S(gJ, thf),
                    'f zeros, phase theta_f': S(gf, thf), 'f-bar zeros, phase theta_f': S(gfb, thf)}
def cj(x):
    if isinstance(x, complex): return [x.real, x.imag]
    if isinstance(x, dict): return {k: cj(v) for k, v in x.items()}
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    return x
out = dict(predictions_sha256=h, per_L=cj(res), pooled=pool, cross=cj(cross))
json.dump(out, open('results.json', 'w'), indent=1)
# table
T = ['| L | family | d | Q | eps | c_L | T2 | zeros | complete (drift, added) | predicted sum | measured sum | slope ± se_block (se_heur) | Hann slope ± se_ctrl (post hoc) | pass (registered rule) |', '|' + '---|' * 14]
for key, r in res.items():
    e = P['Lfunctions'][key]; ci = r['completeness']
    comp = '%+.2f, %d' % (ci['drift'], ci['added_by_gap_search'])
    pr = complex(*e['pred_sharp']); G = complex(*r['G'])
    if 'slope' in r:
        sl = '%.3f ± %.3f (%.3f)' % (r['slope'], r['se_block'], r['se_heur']); hs = '%.3f ± %.3f' % (r['slope_hann'], r['slope_hann_se_ctrl'])
    else:
        sl = 'c_hat %+.3f ± %.3f (%.3f) [naive %.3f]' % (r['c_hat'], r['se_block'], r['se_heur'], r['c_naive']); hs = 'c_hat %+.3f ± %.3f' % (r['c_hat_hann'], r['c_hat_hann_se_ctrl'])
    T.append('| %s | %s | %d | %g | %.3f%+.3fi | %.4f | %d | %d | %s | %.1f%+.1fi | %.1f%+.1fi | %s | %s | %s |' % (
        key, e['family'], e['degree'], e['Q'], e['eps'][0], e['eps'][1], e['c'], e['T2'], r['n'], comp, pr.real, pr.imag, G.real, G.imag, sl, hs,
        'yes' if r['PASS'] else 'NO'))
T.append('')
T.append('| pooled | slope | se_block | se_heur |'); T.append('|---|---|---|---|')
for k, v in pool.items():
    if k == 'lines': continue
    T.append('| %s | %.4f | %.4f | %.4f |' % (k, v['slope'], v['se_block'], v['se_heur']))
T.append(''); T.append('| lines r | pooled slope | se_heur |'); T.append('|---|---|---|')
for k, v in ln.items(): T.append('| %s | %.3f | %.3f |' % (k, v['slope'], v['se_heur']))
T.append('')
for k, v in cross.items():
    T.append('Cross-structure %s: counts %s' % (k, v['counts']))
    T.append(''); T.append('| reading | predicted | measured |'); T.append('|---|---|---|')
    for kk, vv in v.items():
        if kk == 'counts': continue
        pv = P['cross_structure'][k][kk]
        pv = complex(*pv) if isinstance(pv, list) else complex(pv)
        T.append('| %s | %.1f%+.1fi | %.1f%+.1fi |' % (kk, pv.real, pv.imag, vv.real, vv.imag))
    T.append('')
open('results_table.md', 'w').write('\n'.join(T) + '\n')
print('\n'.join(T))
