import os, json, math, hashlib, numpy as np
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
BL = os.path.dirname(os.path.abspath(__file__))+'/../../../../../../../writeup/code/study/lowzero/universal/blind/'
PB = json.load(open(BL + 'predictions_blind.json'))
raw = json.load(open(OUT + 'a1_raw.json')); conv = json.load(open(OUT + 'a2_conv.json')); lev = json.load(open(OUT + 'a3_level.json'))
keys = sorted(raw, key=lambda k: -raw[k]['d'])
diff = []
for k in keys:
    P = {int(A): complex(re, im) for A, re, im in conv[k]['partial']}
    diff.append((P[30] - P[12]).real)
sig_tr = float(np.sqrt(np.mean(np.square(diff))))
fields = {}
reg = []; am = []; sd = []; sdr = []; amr = []
for k in keys:
    r = raw[k]; e = PB['products'][k]
    P = {int(A): [re, im] for A, re, im in conv[k]['partial']}
    sd_reg = e['readings']['M1']['sd']; sd_am = math.sqrt(sd_reg ** 2 + sig_tr ** 2)
    Mr = e['rho']['M']; sd_rho = math.sqrt(PB['noise_kappa']['prod'] * Mr)
    fields[k] = dict(d=r['d'], polynomial=r['pol'], c_rho=r['c_rho_registered'],
                     M1_registered=r['M1_registered'], M1_amended=r['M1_amended'], M1_amended_hann=r['M1_amended_hann'],
                     rho_zeros_theta_K_amended=r['M1_rho_part_amended'], zeta_zeros_theta_K_M3_amended=0.0,
                     M1_amended_by_a=r['by_a'], M1_partial_sums_A_le={str(A): v for A, v in P.items()},
                     sd_M1_registered=sd_reg, sd_truncation=sig_tr, sd_M1_amended=sd_am, sd_rho_part=sd_rho,
                     sd_M1_hann_registered=e['readings']['M1']['sd_hann'])
    reg.append(r['M1_registered']); am.append(r['M1_amended'][0]); sd.append(sd_am); sdr.append(sd_reg)
reg = np.array(reg); am = np.array(am); sd = np.array(sd); sdr = np.array(sdr)
# registered D8 reading if the amended mechanism is true
b_D8 = float((am * reg).sum() / (reg ** 2).sum()); se_D8 = float(math.sqrt((reg ** 2 * sdr ** 2).sum()) / (reg ** 2).sum())
# proposed reading on the amended prediction
se_new = float(math.sqrt((am ** 2 * sd ** 2).sum()) / (am ** 2).sum())
b_new_if_registered_true = float((reg * am).sum() / (am ** 2).sum())
# two-component regression: Re M1 = b1 reg + b2 (am - reg)
X = np.vstack([reg, am - reg]).T; W = 1 / sd ** 2
A = X.T @ (X * W[:, None]); cov = np.linalg.inv(A)
summary = dict(
    T2=300, n_fields=len(keys), truncation='lines r = a/b with a <= 12, all b <= 1.25 a T/2pi + 3, exact integrals (1/2pi) int_0^T e^{2i theta_zeta(t)} r^{it} dt',
    sigma_truncation_rms_S30_minus_S12=sig_tr,
    registered_D8_expected_slope_if_amended_true=b_D8, registered_D8_se=se_D8,
    proposed_D8a_slope_on_amended_expected=1.0, proposed_D8a_se=se_new,
    proposed_D8a_expected_if_registered_mechanism_true=b_new_if_registered_true,
    proposed_D8a_H_prod0_power_sigma=1.0 / se_new,
    regression_on_registered_and_increment_se=[float(math.sqrt(cov[0, 0])), float(math.sqrt(cov[1, 1]))],
    rms_registered=float(np.sqrt(np.mean(reg ** 2))), rms_amended=float(np.sqrt(np.mean(am ** 2))),
    corr_registered_amended=float(np.corrcoef(reg, am)[0, 1]),
    d239=dict(M1_registered=fields['T5_d239']['M1_registered'], M1_amended=fields['T5_d239']['M1_amended'],
              partial=fields['T5_d239']['M1_partial_sums_A_le'], sd=fields['T5_d239']['sd_M1_amended']))
doc = dict(title='Amendment to the round-151 product readings (T5): complete line mechanism for M1', basis_predictions_sha256='c8a77047709d6e8c4fe5dbfb1fa595961dc3ef25ef9d2082ea7245eece3bf50f',
           computed_from='predictions_blind.json constants only (c_rho, eps_rho, field polynomials); no zero of any test L-function read',
           summary=summary, fields=fields,
           level_type=dict(T1_exact_prime_side=lev['T1'], T1_pooled=lev['T1_pooled'], T5_d484_Artin_R=lev['T5_d484_R'],
                           note='No registered L of degree d >= 2 has Q a perfect d-th power (T2 squarefree N >= 3; T3 Q = N^2, d = 3; T4a Q = N^3, d = 4; T4b N in 349..997, no 4th power; T5 zeta_K Q = |d| no cube; T5g Q = 23^2, 31^2 with d = 4). Only L(rho) of d = -484 (degree 2, Q = 22^2) and T1 (degree 1) have coherent prime-power terms.'))
s = json.dumps(doc, indent=1, sort_keys=True)
open(OUT + 'amendment_products.json', 'w').write(s)
h = hashlib.sha256(s.encode()).hexdigest()
open(OUT + 'amendment_products.sha256', 'w').write(h + '  amendment_products.json\n')
print('sha256', h)
print(json.dumps(summary, indent=1))
for k in keys:
    f = fields[k]
    print('%-8s d=%5d  registered %+7.2f  amended %+7.2f%+6.2fi  (A<=3 %+6.1f, A<=30 %+6.1f)  sd %.1f' % (k, f['d'], f['M1_registered'], f['M1_amended'][0], f['M1_amended'][1], f['M1_partial_sums_A_le']['3'][0], f['M1_partial_sums_A_le']['30'][0], f['sd_M1_amended']))
# T1 complex pooled of exact over law
num = sum(complex(*v['exact_prime_side']) * complex(*v['pred_registered']).conjugate() for v in lev['T1'].values())
den = sum(abs(complex(*v['pred_registered'])) ** 2 for v in lev['T1'].values())
print('T1 complex pooled exact/law: %.4f%+.4fi' % ((num / den).real, (num / den).imag))
dev = [abs(complex(*v['exact_prime_side']) - complex(*v['pred_registered'])) / v['sd_along'] for v in lev['T1'].values()]
print('T1 |exact - law| / sd_along: median %.2f, max %.2f, n > 3: %d' % (np.median(dev), max(dev), sum(1 for x in dev if x > 3)))
