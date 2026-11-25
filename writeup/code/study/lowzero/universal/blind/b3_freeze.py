"""freeze the blind predictions: reads <testset.json> and <constants dir>, writes <out.json> and prints its sha256.
Noise model (registered; measured in round 1 on 20 primitive L, sharp sums, Var/M): r = 1 along the prediction 0.077,
across 0.21 (Hann: the same factors per 3M/8); lines along 0.14, across 0.18; products at a non-own frequency (no rigidity):
1.0; rho zeros under the zeta phase (form factor tau ~ 0.4): 0.4. Expected errors below follow from these.
usage: b3_freeze.py testset.json constants out.json"""
import sys, os, json, math, cmath, hashlib, datetime, glob
import numpy as np
TS = json.load(open(sys.argv[1])); CD = sys.argv[2]; OUT = sys.argv[3]
K = dict(r1_along=0.077, r1_across=0.21, line_along=0.14, line_across=0.18, prod=1.0, zeta_phase_on_rho=0.4)
H = TS['height']
C = {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(os.path.join(CD, '*.json')))}
ORDER = ['T5', 'T5g', 'T1', 'T2', 'T3', 'T4a', 'T4b']
EST = dict(T1=45, T2=30, T3=60, T4a=60, T4b=60, T5=10, T5g=60)       # seconds per L (proxy timings, b1_full_out_*.txt)
P = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), seed=TS['seed'],
         testset_sha256=hashlib.sha256(open(sys.argv[1], 'rb').read()).hexdigest(), heights=H, noise_kappa=K,
         law='sum_{0<g<=T2} e^{2i theta_L(g)} r^{-ig} = -eps c_L rho(r)/sqrt(ab) T2/2pi; c_L = B_L(1)/L(1, Ad); products c = 0, '
             'with the finite-height mechanism M1 = -c_rho sum_{k<=T2/2pi} rho(1/k), M2 = -psi_rho(T2/2pi), M3 = 0',
         manifest=[dict(key='Zeta', cls='aux', spec=['zeta'], T2=300, poles=1, est_secs=5)], L={}, products={}, pairs=[])
def sd(kappa, M): return math.sqrt(kappa * M)
for cls in ORDER:
    for key, c in C.items():
        if c.get('cls') != cls: continue
        T2 = c['T2']; unit = T2 / (2 * math.pi); eps = complex(*c['eps']); M = c['M']
        P['manifest'].append(dict(key=key, cls=cls, spec=c['spec'], T2=T2, poles=1 if c['spec'][0] == 'nf' else 0, est_secs=EST[cls]))
        if cls in ('T5', 'T5g'):
            e = dict(cls=cls, spec=c['spec'], Q=c['Q'], mu=c['mu'], eps=c['eps'], T2=T2, M=M, c=0.0)
            if cls == 'T5':
                r = c['rho']; m = c['mech']; Mr = r['M']; Mz = 138.0      # zeros of zeta below 300
                e.update(d=c['d'], rho=r, mech=m,
                         readings=dict(
                             M1=dict(pred=m['M1'], pred_hann=m['M1_hann'], H_prod0=0.0, sd=sd(K['prod'], M), sd_hann=sd(K['prod'], 3 * M / 8)),
                             M2=dict(pred=m['M2'], pred_hann=m['M2_hann'], H_prod0=0.0, sd=sd(K['zeta_phase_on_rho'], Mr), sd_hann=sd(K['zeta_phase_on_rho'], 3 * Mr / 8)),
                             M3=dict(pred=0.0, sd=sd(K['prod'], Mz)),
                             R_artin=dict(pred=-complex(*r['eps']).real * r['c'] * unit, pred_hann=-complex(*r['eps']).real * r['c'] * unit / 2,
                                          H0=0.0, H_c1=-unit, sd=sd(K['r1_along'], Mr), sd_hann=sd(K['r1_along'], 3 * Mr / 8),
                                          lines={kk: dict(pred=[v[0] * unit, v[1] * unit], sd=sd(K['line_along'], Mr)) for kk, v in r['lines'].items()})))
            else:
                hp = complex(*c['half_phase_pred'])
                e.update(factors=c['factors'], readings=dict(half_phase=dict(pred=hp.real, pred_hann=hp.real / 2, H0=0.0, sd=sd(K['r1_along'], M)),
                                                             product_phase=dict(pred=0.0, sd=sd(K['prod'], M))))
            P['products'][key] = e
            continue
        cc = c['c']; G = -eps * cc * unit
        e = dict(cls=cls, spec=c['spec'], Q=c['Q'], mu=c['mu'], degree=c['degree'], eps=c['eps'], c=cc, c_method=c['c_method'],
                 c_direct=c.get('c_direct'), T2=T2, M=M, pred=[G.real, G.imag], pred_hann=[G.real / 2, G.imag / 2],
                 alt=dict(H0=[0, 0], H_noroot=[-cc * unit, 0], H_conj=[(-eps.conjugate() * cc * unit).real, (-eps.conjugate() * cc * unit).imag],
                          H_c1=[(-eps * unit).real, (-eps * unit).imag]),
                 sd=dict(along=sd(K['r1_along'], M), across=sd(K['r1_across'], M), along_hann=sd(K['r1_along'], 3 * M / 8),
                         line_along=sd(K['line_along'], M), line_across=sd(K['line_across'], M)),
                 lines={kk: dict(pred=[v[0] * unit, v[1] * unit], pred_hann=[v[0] * unit / 2, v[1] * unit / 2]) for kk, v in c['lines'].items()})
        P['L'][key] = e
for a, b in TS['T5pairs']['test']:
    def nw(key):
        parts = key.split('_'); N = int(parts[0][1:]); return N, [int(x.replace('m', '-')) for x in parts[1:]]
    P['pairs'].append(dict(keys=[a, b], N=[nw(a)[0], nw(b)[0]], ainv=[nw(a)[1], nw(b)[1]], T2=H['T5pairs'], pred=0.0, pred_hann=0.0,
                           note='union of the two zero lists (round 144), phase theta_1 + theta_2; sd = sqrt(n1 + n2) (no rigidity)'))
# pooled expected errors (registered)
def pooled(keys):
    num = sum(abs(complex(*P['L'][k]['pred'])) ** 2 for k in keys)
    var = sum(abs(complex(*P['L'][k]['pred'])) ** 2 * P['L'][k]['sd']['along'] ** 2 for k in keys)
    return math.sqrt(var) / num if num > 0 else None
fam = {}
for k, e in P['L'].items(): fam.setdefault(e['cls'], []).append(k)
P['expected_se'] = {f: pooled(v) for f, v in fam.items()}
P['expected_se']['ALL'] = pooled(list(P['L']))
lnum = sum(abs(complex(*l['pred'])) ** 2 for e in P['L'].values() for l in e['lines'].values())
lvar = sum(abs(complex(*l['pred'])) ** 2 * e['sd']['line_along'] ** 2 for e in P['L'].values() for l in e['lines'].values())
P['expected_se']['lines'] = math.sqrt(lvar) / lnum
for rd in ('M1', 'M2'):
    num = sum(e['readings'][rd]['pred'] ** 2 for e in P['products'].values() if e['cls'] == 'T5')
    var = sum(e['readings'][rd]['pred'] ** 2 * e['readings'][rd]['sd'] ** 2 for e in P['products'].values() if e['cls'] == 'T5')
    P['expected_se']['T5_' + rd] = math.sqrt(var) / num if num > 0 else None
num = sum(e['readings']['R_artin']['pred'] ** 2 for e in P['products'].values() if e['cls'] == 'T5')
var = sum(e['readings']['R_artin']['pred'] ** 2 * e['readings']['R_artin']['sd'] ** 2 for e in P['products'].values() if e['cls'] == 'T5')
P['expected_se']['T5_R_artin'] = math.sqrt(var) / num if num > 0 else None
P['decisions'] = {
 'D1 family': 'pooled slope b_F = sum Re(G conj pred)/sum |pred|^2 over the COMPLETE L of the family; PASS iff |b_F - 1| <= 3 se_F (se_F from the noise model, frozen in expected_se, recomputed on the complete subset)',
 'D2 all primitive': 'same over all primitive L; PASS iff |b - 1| <= 3 se',
 'D3 H0': 'H0 (no main term) rejected iff b/se >= 5 (all primitive)',
 'D4 lines': 'pooled slope over r = 2, 1/2, 3, 1/3 and all primitive L (line noise); PASS iff |b - 1| <= 3 se',
 'D5 amplitude follows c': 'weighted LS of c_hat_L = Re(G conj(-eps))/(T2/2pi) on [1, c_L - 1] over all primitive L; PASS iff |beta - 1| <= 3 se_beta; H-c1 (beta = 0) rejected iff beta/se_beta >= 3',
 'D6 root number eps = -1': 'pooled slope over the primitive L with eps = -1; PASS iff |b - 1| <= 3 se (H-noroot: -1)',
 'D7 complex eps (T1)': 'complex pooled slope s = sum G conj pred/sum |pred|^2; PASS iff |Im s| <= 3 se_im (H-conj: Im s = registered value)',
 'D8 products M1': 'b = sum M1 mech/sum mech^2 over the 30 fields; PASS (mechanism) iff |b - 1| <= 3 se; H-prod0 rejected iff b/se >= 3',
 'D9 products M2': 'same with -psi_rho; M3: |mean| <= 3 sd/sqrt(n)',
 'D10 Artin law on the rho zeros of T5': 'pooled slope; PASS iff |b - 1| <= 3 se',
 'D11 GL2-type genus-2 products (T5g)': 'half-phase slope within 3 se of 1 (pooled over the 2 curves); product phase |c_hat| <= 3 se',
 'D12 pairs': 'mean of the 10 product-phase sums within 3 sd/sqrt(10) of 0',
 'falsifier': 'the general law is falsified if the slope pooled over all primitive test L lies outside 1 +- 4 se (se = expected_se ALL), or any family slope outside 1 +- 4 se_F'}
s = json.dumps(P, indent=1, sort_keys=True)
open(OUT, 'w').write(s)
print('L:', len(P['L']), 'products:', len(P['products']), 'pairs:', len(P['pairs']), 'manifest:', len(P['manifest']))
print('expected se:', {k: (round(v, 4) if v else v) for k, v in P['expected_se'].items()})
print('sha256(%s) = %s' % (OUT, hashlib.sha256(s.encode()).hexdigest()))
