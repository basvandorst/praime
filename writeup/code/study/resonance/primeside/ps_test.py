"""ps_test.py: readings and registered decisions of the prime-side test.  Refuses any predictions file but the frozen one.

usage: ps_test.py RUN        RUN = dry (seen range, not blind) or fresh (registered)
Reads predictions_primeside.json (sha256 must equal FROZEN) and sums_RUN.jsonl (ps_sums.py), writes readings_RUN.json."""
import sys, os, json, math, hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FROZEN = '63e00dc97e51d856504b15f50b7d96ae6d63819d7082e8c1b3e3dc414e56c790'

def C(v): return complex(v[0], v[1])

def aggregate(run, pred):
    """sum the chunk records over the readings range of each series; cumulative sharp sums on the grid."""
    recs = [json.loads(l) for l in open(os.path.join(HERE, 'sums_%s.jsonl' % run))]
    agg, zones, curve = {}, {}, {}
    for key, ent in pred['series'].items():
        x0, x1 = ent['runs'][run]['read']; grid = ent['runs'][run]['grid']
        sharp = np.zeros(9, complex); smooth = np.zeros((9, 9), complex); nsh = 0.0; nsm = np.zeros(9); cnt = 0
        cum = {g: 0j for g in grid}; seen = set()
        for r in recs:
            if key not in r['res'] or r['L'] < x0 or r['U'] > x1: continue
            assert (r['L'], r['U']) not in seen; seen.add((r['L'], r['U']))
            q = r['res'][key]
            sharp += np.array([C(z) for z in q['sharp']]); smooth += np.array([[C(z) for z in row] for row in q['smooth']])
            nsh += q['noise_sharp']; nsm += np.array(q['noise_smooth']); cnt += q['count']
            for g in grid:
                if r['U'] <= g: cum[g] += C(q['sharp'][0])
        covered = sum(b - a for (a, b) in seen)
        agg[key] = dict(sharp=sharp, smooth=smooth, noise_sharp=nsh, noise_smooth=nsm, count=cnt, covered=covered, span=x1 - x0, cum=[cum[g] for g in grid])
    x0, x1 = pred['series']['A4']['runs'][run]['read']
    for r in recs:
        if r['group'] != 'ABC' or r['L'] < x0 or r['U'] > x1 or r.get('zones') is None: continue
        for name, v in r['zones'].items():
            z = zones.setdefault(name, {'0': 0j, '4': 0j, 'noise': 0.0})
            z['0'] += C(v['0']); z['4'] += C(v['4']); z['noise'] += v['noise']
    return agg, zones

def readings(run, pred):
    agg, zones = aggregate(run, pred); dec = pred['decisions']; out = {'run': run, 'series': {}, 'zones': {}}
    for key, ent in pred['series'].items():
        a = agg[key]; pr = ent['runs'][run]; c = ent['c']
        if a['covered'] != a['span']:
            out['series'][key] = {'reading': 'no reading (range not covered: %d of %d)' % (a['covered'], a['span'])}; continue
        MTs, PPs = C(pr['MT_sharp']), C(pr['PP_sharp'])
        MTw = [C(z) for z in pr['MT_smooth']]; PPw = [C(z) for z in pr['PP_smooth']]
        rs = (a['sharp'][0] + PPs) / MTs
        ss = math.sqrt(np.mean(np.abs(a['sharp'][1:]) ** 2) / 2) / abs(MTs)
        rw = [(a['smooth'][i][0] + PPw[i]) / MTw[i] for i in range(9)]
        sw = [math.sqrt(np.mean(np.abs(a['smooth'][i][1:]) ** 2) / 2) / abs(MTw[i]) for i in range(9)]
        rand_s = math.sqrt(a['noise_sharp'] / 2) / abs(MTs); rand_w = math.sqrt(a['noise_smooth'][0] / 2) / abs(MTw[0])
        sub = np.array([z.real for z in rw[1:]]); subs = np.array(sw[1:])
        wm = np.sum(sub / subs ** 2) / np.sum(1 / subs ** 2)
        chi_sub = float(np.sum(((sub - wm) / subs) ** 2) / (len(sub) - 1))
        noPP_w = a['smooth'][0][0] / MTw[0]
        out['series'][key] = dict(c=c, K=ent['K'], count=a['count'],
            smooth=[rw[0].real, rw[0].imag], smooth_err=sw[0], smooth_err_random_model=rand_w,
            sharp=[rs.real, rs.imag], sharp_err=ss, sharp_err_random_model=rand_s,
            sub_windows=[[z.real, z.imag] for z in rw[1:]], sub_errs=sw[1:], sub_mean=wm, sub_chi2_dof=chi_sub,
            ctrl_smooth=[[z.real / abs(MTw[0]), z.imag / abs(MTw[0])] for z in a['smooth'][0][1:]],
            smooth_without_PP=[noPP_w.real, noPP_w.imag], PP_over_MT=abs(PPw[0] / MTw[0]),
            alt_constant1=1 / c,
            sharp_curve=[[g, (z + C(pp)).real / abs(C(mt)) if abs(C(mt)) > 0 else 0, (z + C(pp)).imag / abs(C(mt)) if abs(C(mt)) > 0 else 0]
                         for g, z, mt, pp in zip(pr['grid'], a['cum'], pr['MT_grid'], pr['PP_grid'])])
    # pooled readings (smooth, registered)
    keys = [k for k, v in out['series'].items() if 'smooth' in v]
    def pool(ks, which='smooth'):
        x = np.array([out['series'][k][which][0] for k in ks]); s = np.array([out['series'][k][which + '_err'] for k in ks])
        w = 1 / s ** 2; m = float(np.sum(w * x) / np.sum(w)); e = float(1 / math.sqrt(np.sum(w)))
        chi = float(np.sum(((x - m) / s) ** 2) / max(1, len(ks) - 1))
        return dict(n=len(ks), mean=m, err=e, chi2_dof=chi, err_scaled=e * math.sqrt(max(1.0, chi)))
    fam = {'all': keys, 'QI_A': [k for k in keys if k[0] == 'A'], 'QS3_B': [k for k in keys if k[0] == 'B'],
           'dihedral_C': [k for k in keys if k[0] == 'C'], 'non_CM_D_E': [k for k in keys if k[0] in 'DE'], 'degree2': [k for k in keys if k != 'E11a']}
    out['pooled'] = {f: pool(ks) for f, ks in fam.items() if ks}
    out['pooled_sharp'] = {f: pool(ks, 'sharp') for f, ks in fam.items() if ks}
    # alternatives: chi^2 over all series (smooth, both components)
    chi = {'law': 0.0, 'constant_1': 0.0, 'no_effect': 0.0, 'law_without_PP': 0.0}
    for k in keys:
        v = out['series'][k]; s = v['smooth_err']; x = complex(*v['smooth'])
        chi['law'] += abs(x - 1) ** 2 / s ** 2; chi['constant_1'] += abs(x - v['alt_constant1']) ** 2 / s ** 2
        chi['no_effect'] += abs(x) ** 2 / s ** 2; chi['law_without_PP'] += abs(complex(*v['smooth_without_PP']) - 1) ** 2 / s ** 2
    out['chi2'] = chi; out['chi2_n'] = 2 * len(keys)
    # controls: all control readings in units of the main term (should scatter with sd = smooth_err)
    cz = np.array([z for k in keys for z in out['series'][k]['ctrl_smooth']]) / np.repeat([out['series'][k]['smooth_err'] for k in keys], 8)[:, None]
    out['controls'] = dict(n=int(len(cz)), rms_in_errors=float(np.sqrt(np.mean(cz ** 2))))
    # zones
    zp = pred['zones'][run]
    for name, z in zones.items():
        p0, p4 = C(zp[name]['0']), C(zp[name]['4']); sd = math.sqrt(z['noise'] / 2)
        out['zones'][name] = dict(m0=[z['0'].real, z['0'].imag], pred_m0=[p0.real, p0.imag], m4=[z['4'].real, z['4'].imag], pred_m4=[p4.real, p4.imag],
                                  sd=sd, dev_m0=abs(z['0'] - p0) / sd, dev_m4=abs(z['4'] - p4) / sd,
                                  ratio_m0=[(z['0'] / p0).real, (z['0'] / p0).imag] if abs(p0) > 5 * sd else None)
    # Pythagorean pooled amplitude: projection of measured (m = 0) on predicted over the Pythagorean windows with nonzero weight
    pyth = [n for n in zones if n[0].isdigit() and abs(C(zp[n]['0'])) > 3 * math.sqrt(zones[n]['noise'] / 2)]
    if pyth:
        num = sum((zones[n]['0'] * C(zp[n]['0']).conjugate()).real / (zones[n]['noise'] / 2) for n in pyth)
        den = sum(abs(C(zp[n]['0'])) ** 2 / (zones[n]['noise'] / 2) for n in pyth)
        out['zones_pooled_pyth'] = dict(windows=pyth, amplitude=num / den, err=1 / math.sqrt(den))
    return out

def decide(out, dec):
    """the registered decisions (thresholds from decisions.json)."""
    D = {}
    P = out['pooled']['all']; e = max(P['err'], P['err_scaled'])
    D['pooled slope within 1 +- %.2f and within %g errors' % (dec['pooled_band'], dec['pooled_nsig'])] = abs(P['mean'] - 1) <= dec['pooled_band'] and abs(P['mean'] - 1) <= dec['pooled_nsig'] * e
    for f in ('QI_A', 'QS3_B', 'dihedral_C', 'non_CM_D_E'):
        if f in out['pooled']:
            Q = out['pooled'][f]; e = max(Q['err'], Q['err_scaled'])
            D['family %s within %g errors of 1' % (f, dec['family_nsig'])] = abs(Q['mean'] - 1) <= dec['family_nsig'] * e
    nbad = sum(1 for k, v in out['series'].items() if 'smooth' in v and abs(complex(*v['smooth']) - 1) > dec['series_nsig'] * v['smooth_err'])
    D['at most %d series with |reading - 1| > %g errors' % (dec['series_max_bad'], dec['series_nsig'])] = nbad <= dec['series_max_bad']
    ch = out['chi2']
    D['law beats constant 1 by Delta chi2 >= %g' % dec['delta_chi2']] = ch['constant_1'] - ch['law'] >= dec['delta_chi2']
    D['law beats no effect by Delta chi2 >= %g' % dec['delta_chi2']] = ch['no_effect'] - ch['law'] >= dec['delta_chi2']
    D['controls rms in errors within [%g, %g]' % tuple(dec['controls_rms'])] = dec['controls_rms'][0] <= out['controls']['rms_in_errors'] <= dec['controls_rms'][1]
    if 'zones_pooled_pyth' in out:
        Z = out['zones_pooled_pyth']
        D['Pythagorean rays: pooled amplitude within 1 +- %g errors' % dec['zone_nsig']] = abs(Z['amplitude'] - 1) <= dec['zone_nsig'] * Z['err']
        zero = [n for n in dec['zone_zero_windows'] if n in out['zones']]
        D['zero-weight and control windows: all within %g sd of the ray model' % dec['zone_nsig']] = all(out['zones'][n]['dev_m0'] <= dec['zone_nsig'] for n in zero)
    falsified = (abs(out['pooled']['all']['mean'] - 1) > dec['falsifier_band']) or (ch['law'] > ch['constant_1'])
    return D, falsified

def report(out, D, falsified):
    L = []
    L.append('RUN %s' % out['run'])
    L.append('%-5s %9s %22s %7s %22s %7s %7s %9s %9s' % ('case', 'c', 'smooth reading', 'err', 'sharp reading', 'err', 'rand.m', 'sub chi2', 'const 1'))
    for k, v in out['series'].items():
        if 'smooth' not in v: L.append('%-5s %s' % (k, v['reading'])); continue
        L.append('%-5s %9.5f %10.4f %+10.4fi %7.4f %10.4f %+10.4fi %7.4f %7.4f %9.2f %9.3f' % (k, v['c'], v['smooth'][0], v['smooth'][1], v['smooth_err'],
                 v['sharp'][0], v['sharp'][1], v['sharp_err'], v['smooth_err_random_model'], v['sub_chi2_dof'], v['alt_constant1']))
    for f, P in out['pooled'].items():
        L.append('pooled %-10s n=%2d  slope %.4f +- %.4f (chi2/dof %.2f; scaled err %.4f)   sharp: %.4f +- %.4f' % (f, P['n'], P['mean'], P['err'], P['chi2_dof'], P['err_scaled'],
                 out['pooled_sharp'][f]['mean'], out['pooled_sharp'][f]['err']))
    L.append('chi2 (%d readings): ' % out['chi2_n'] + '  '.join('%s %.1f' % kv for kv in out['chi2'].items()))
    L.append('controls: %d readings, rms %.3f errors' % (out['controls']['n'], out['controls']['rms_in_errors']))
    for n, z in out['zones'].items():
        L.append('zone %-7s m=0: %10.0f %+10.0fi  model %10.0f %+10.0fi  (%.1f sd)   m=4: %10.0f %+10.0fi  model %10.0f %+10.0fi (%.1f sd)  sd %.0f' % (
            n, z['m0'][0], z['m0'][1], z['pred_m0'][0], z['pred_m0'][1], z['dev_m0'], z['m4'][0], z['m4'][1], z['pred_m4'][0], z['pred_m4'][1], z['dev_m4'], z['sd']))
    if 'zones_pooled_pyth' in out:
        Z = out['zones_pooled_pyth']; L.append('Pythagorean windows %s: pooled amplitude %.4f +- %.4f' % (','.join(Z['windows']), Z['amplitude'], Z['err']))
    for d, v in D.items(): L.append(('PASS ' if v else 'FAIL ') + d)
    L.append('FALSIFIER MET' if falsified else 'falsifier not met')
    return '\n'.join(L)

def main():
    run = sys.argv[1]
    fn = os.path.join(HERE, 'predictions_primeside.json'); s = open(fn).read()
    h = hashlib.sha256(s.encode()).hexdigest()
    if h != FROZEN: sys.exit('refused: predictions_primeside.json has sha256 %s, frozen %s' % (h, FROZEN))
    pred = json.loads(s); out = readings(run, pred); D, fals = decide(out, pred['decisions'])
    out['decisions'] = {k: bool(v) for k, v in D.items()}; out['falsified'] = bool(fals)
    def enc(o):
        if isinstance(o, complex): return [o.real, o.imag]
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.floating,)): return float(o)
        raise TypeError(type(o))
    json.dump(out, open(os.path.join(HERE, 'readings_%s.json' % run), 'w'), indent=1, default=enc)
    print(report(out, D, fals))

if __name__ == '__main__':
    main()
