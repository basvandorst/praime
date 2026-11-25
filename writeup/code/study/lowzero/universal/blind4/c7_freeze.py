"""freeze the degree-4 predictions. usage: c7_freeze.py <constants dir> <out.json> [<tag>]
Per L: sharp sum on (0, T2], Hann (sin^2(pi g/T2)) half of it, band (T2/4, T2]: -eps c (3 T2/4)/2pi, lines r = 2, 1/2, 3, 1/3;
alternatives H0 (0), H_noroot (-c T2/2pi), H_c1 (-eps T2/2pi), H_eps0 (Sym^3 with eps = -1 carry no main term: 0).
Noise for EXPECTED errors only (measured in round 151, Var/M): Sym^3 along 0.479, across 0.504, Hann 0.464 per 3M/8, lines 0.353/0.350;
genus 2 0.152, 0.175, Hann 0.295, lines 0.214/0.324; Rankin-Selberg (not measured): the Sym^3 values. Decisions use jackknife errors.
Cost estimate per L (seconds, timings of this folder incl. a 1.3 factor for the gap re-search): Sym^3 S70 1.3*max(20, 126 (N/142)^1.5),
S50 1.3*max(20, 76 (N/215)^1.5); genus 2 1.3*(15 + 90 sqrt(N/2981)); RS 1.3*12.5*(N1 N2/266) at T2 = 60 (measured 12.5 s at 266, 69 s at 1443)."""
import sys, os, json, math, glob, hashlib, datetime
CD, OUT = sys.argv[1], sys.argv[2]
K = dict(sym3=dict(along=0.479, across=0.504, hann=0.464, l_along=0.353, l_across=0.350),
         g2=dict(along=0.152, across=0.175, hann=0.295, l_along=0.214, l_across=0.324),
         rs=dict(along=0.479, across=0.504, hann=0.464, l_along=0.353, l_across=0.350))
def est(c):
    s = c['spec']
    if s[0] == 'sym':
        N = c['Q'] ** (1 / 3)
        return 1.3 * (max(20, 126 * (N / 142) ** 1.5) if c['T2'] == 70 else max(20, 76 * (N / 215) ** 1.5))
    if s[0] == 'g2': return 1.3 * (15 + 90 * math.sqrt(c['Q'] / 2981))
    return 1.3 * 12.5 * (c['T2'] / 60) ** 3 * math.sqrt(c['Q']) / 266
P = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), noise_kappa=K,
         law='sum_{0<g<=T2} e^{2i theta_L(g)} r^{-ig} = -eps c_L rho_L(r)/sqrt(ab) T2/2pi; c_L = B_L(1)/L(1, Ad L)',
         completeness='bz4_lib.py (corrected offset c_exp = 1/2 + poles - m0/2 - A0/pi)', L={}, manifest=[])
for f in sorted(glob.glob(os.path.join(CD, '*.json'))):
    c = json.load(open(f)); T2 = c['T2']; unit = T2 / (2 * math.pi); eps = complex(*c['eps']); cc = c['c']; M = c['M']
    cls = c['cls']; k = K[cls]; G = -eps * cc * unit; Gb = -eps * cc * (0.75 * T2) / (2 * math.pi)
    e = dict(cls=cls, stratum=('S70' if T2 == 70 else 'S50') if cls == 'sym3' else None, spec=c['spec'], Q=c['Q'], mu=c['mu'],
             eps=c['eps'], c=cc, c_direct_end=c['c_direct_end'], c_direct_P3=c['c_direct_P3'], c_exact_check=c.get('c_exact_check'),
             T2=T2, M=M, M_band=c['M_band'], pred=[G.real, G.imag], pred_hann=[G.real / 2, G.imag / 2], pred_band=[Gb.real, Gb.imag],
             lines={nm: [v[0] * unit, v[1] * unit] for nm, v in c['lines'].items()},
             alt=dict(H0=0.0, H_noroot=-cc * unit, H_c1=(-eps * unit).real, H_eps0=0.0 if (cls == 'sym3' and eps.real < 0) else G.real),
             sd=dict(along=math.sqrt(k['along'] * M), across=math.sqrt(k['across'] * M), hann=math.sqrt(k['hann'] * 3 * M / 8),
                     band=math.sqrt(k['along'] * c['M_band']), line_along=math.sqrt(k['l_along'] * M)))
    P['L'][c['key']] = e
    P['manifest'].append(dict(key=c['key'], cls=cls, spec=c['spec'], T2=T2, poles=0, est_secs=round(est(c))))
order = dict(rs=0, g2=1, sym3=2)
P['manifest'].sort(key=lambda m: (order[m['cls']], m['est_secs']))
def model_se(keys, which='along', field='pred'):
    num = sum(abs(complex(*P['L'][k][field])) ** 2 for k in keys)
    var = sum(abs(complex(*P['L'][k][field])) ** 2 * P['L'][k]['sd'][which] ** 2 for k in keys)
    return math.sqrt(var) / num if num else None
sets = {'all': list(P['L']), 'sym3': [k for k in P['L'] if P['L'][k]['cls'] == 'sym3'],
        'sym3_eps-1': [k for k in P['L'] if P['L'][k]['cls'] == 'sym3' and P['L'][k]['eps'][0] < 0],
        'sym3_eps+1': [k for k in P['L'] if P['L'][k]['cls'] == 'sym3' and P['L'][k]['eps'][0] > 0],
        'g2': [k for k in P['L'] if P['L'][k]['cls'] == 'g2'], 'rs': [k for k in P['L'] if P['L'][k]['cls'] == 'rs']}
for s in ('S70', 'S50'):
    for sg, nm in ((-1, '-1'), (1, '+1')):
        sets['sym3_%s_eps%s' % (s, nm)] = [k for k in sets['sym3'] if P['L'][k]['stratum'] == s and P['L'][k]['eps'][0] * sg > 0]
P['sets'] = sets
P['expected_se'] = {nm: model_se(ks) for nm, ks in sets.items() if ks}
P['expected_se']['band_sym3_eps-1'] = model_se(sets['sym3_eps-1'], 'band', 'pred_band')
P['expected_se']['hann_all'] = model_se(sets['all'], 'hann', 'pred_hann')
P['est_total_secs'] = sum(m['est_secs'] for m in P['manifest'])
P['decisions'] = {
 'R1 all degree 4': 'slope b = sum Re(G conj pred)/sum |pred|^2 over the COMPLETE new L; jackknife se over L; PASS iff |b - 1| <= 3 se',
 'R2 Sym^3 by root number (watch item)': 'b(eps = -1), b(eps = +1), and their difference, jackknife se; FLUCTUATION (the law holds for eps = -1) iff |b(-1) - 1| <= 3 se and b(-1)/se >= 3; '
                                          'EFFECT iff 1 - b(-1) > 4 se while |b(+1) - 1| <= 3 se; else UNDECIDED',
 'R2s strata': 'the same four slopes per stratum S70 (T2 = 70, N <= 150) and S50 (T2 = 50, 150 < N <= 230): reported, no decision',
 'R3 genus 2': 'slope, jackknife; PASS iff |b - 1| <= 3 se', 'R4 Rankin-Selberg': 'slope, jackknife; PASS iff |b - 1| <= 3 se',
 'R5 amplitude follows c': 'c_hat_L = Re(G conj(-eps))/(T2/2pi) on [1, c_L - 1] (weights 1/(kappa M)), jackknife over L; PASS iff |beta - 1| <= 3 se; H-c1 rejected iff beta/se >= 3',
 'R6 lines': 'pooled slope over r = 2, 1/2, 3, 1/3 and all L, jackknife over L; PASS iff |b - 1| <= 3 se',
 'R7 Hann': 'pooled Hann slope, jackknife; PASS iff |b - 1| <= 3 se',
 'R8 band (T2/4, T2]': 'slopes for all, Sym^3 eps = -1 and eps = +1 (jackknife); reported; if R2 reads EFFECT and the eps = -1 band slope is within 3 se of 1, the deficit is attributed to the low-height part of the window',
 'R9 with round 151': 'the 30 degree-4 L of round 151 added (their zero files): pooled slopes, reported only',
 'falsifier': 'the law is falsified in degree 4 if the slope pooled over all complete new degree-4 L lies outside 1 +- max(4 se_jk, 0.10)'}
s = json.dumps(P, indent=1, sort_keys=True); open(OUT, 'w').write(s)
print('L', len(P['L']), {k: len(v) for k, v in sets.items()})
print('expected se (measured noise):', {k: round(v, 4) for k, v in P['expected_se'].items() if v})
print('estimated zero cost: %.0f s = %.2f h' % (P['est_total_secs'], P['est_total_secs'] / 3600))
print('sha256(%s) = %s' % (OUT, hashlib.sha256(s.encode()).hexdigest()))
