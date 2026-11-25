"""Freeze per-L-function predictions of Var(s) and <r~> for the reserved zero sets (and the used set R120, for the dry
run).  Only the discriminant lists (the npz key names) and the height bands are read from the reserved files; no zero.
Noise model from the USED set R120 only (residual of data - law after the indicator regression), scaled by sqrt(n).
Output: predictions_spacing.json (+ its sha256 in predictions_spacing.sha256)."""
import sys, os, math, json, hashlib, zipfile
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/lib'); sys.path.insert(0, H)
import numpy as np
from scipy.special import loggamma
import predlib
LZ = os.path.join(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..')), 'data', 'riemann', 'lzeros') + '/'
SETS = {'B': (['r134_B.npz'], 100, 600), 'Dn': (['r136_Dn.npz'], 100, 600), 'Fn': (['r138_Fn.npz'], 100, 600),
        'C': (['r134_C.npz'], 600, 1000), 'Bp': (['r138_Bp.npz'], 600, 1000),
        'R120': (['r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'], 100, 600)}
PR = (2, 3, 5, 7, 11, 13); VARIANTS = ('law', 'law2', 'slow', 'bblm', 'ignored')

def keys_only(f):
    with zipfile.ZipFile(LZ + f) as z:                       # directory of the archive only; no array is read
        return [n[:-4] for n in z.namelist()]

def theta(t, q, a): return np.imag(loggamma(0.25 + a / 2 + 0.5j * t)) + t / 2 * math.log(q / math.pi)

def nexp(D, lo, hi):
    q = abs(D); a = 0 if D > 0 else 1
    return float((theta(hi, q, a) - theta(lo, q, a)) / math.pi)

P = predlib.Predictor(); out = {'note': __doc__, 'primes': PR, 'variants': VARIANTS, 'V0_sine': P.V0, 'R0_sine': P.R0}
for name, (files, lo, hi) in SETS.items():
    Ds = []
    for f in files: Ds += [int(k[1:]) for k in keys_only(f)]
    Ds = sorted(set(Ds), key=lambda d: (abs(d), d))
    rows = {v: {'var': [], 'r': []} for v in VARIANTS}; meanL = []; nsp = []
    for D in Ds:
        pr = P.predict(D, lo, hi, law2=True)
        for v in VARIANTS:
            rows[v]['var'].append(pr[v]['var']); rows[v]['r'].append(pr[v]['r'])
        meanL.append(pr['meanL']); nsp.append(nexp(D, lo, hi) - 1)
    X = np.array([[1.0] + [float(abs(D) % p == 0) for p in PR] for D in Ds])
    coef = {v: {st: np.linalg.lstsq(X, np.array(rows[v][st]), rcond=None)[0].tolist() for st in ('var', 'r')} for v in VARIANTS}
    out[name] = dict(D=Ds, height=[lo, hi], meanL=meanL, n_spacings_expected=nsp, pred=rows, coef_matched=coef,
                     group_sizes=[int(X[:, j].sum()) for j in range(1, 7)], n=len(Ds))
    print('%-4s n=%d heights %d-%d mean L %.3f  spacings/set %.0f  groups %s' % (name, len(Ds), lo, hi, np.mean(meanL), np.mean(nsp), out[name]['group_sizes']))
    for v in VARIANTS:
        print('     %-8s Var level %.5f coef(2,3,5,7) %s | r level %.5f coef %s' % (v, np.mean(rows[v]['var']), ' '.join('%+.5f' % c for c in coef[v]['var'][1:5]), np.mean(rows[v]['r']), ' '.join('%+.5f' % c for c in coef[v]['r'][1:5])))
# noise model from the used set R120: per-L-function residual of (data - law) after the indicator regression
sys.path.insert(0, H); import used_check
Z = used_check.load_R120(); dv, dr, nn = [], [], []
R = out['R120']
for i, D in enumerate(R['D']):
    g = Z['D%d' % D]; g = g[(g >= 100) & (g <= 600)]; v, r, n = used_check.stats(used_check.unfold(g, abs(D)))
    dv.append(v - R['pred']['law']['var'][i]); dr.append(r - R['pred']['law']['r'][i]); nn.append(n)
X = np.array([[1.0] + [float(abs(D) % p == 0) for p in PR] for D in R['D']]); dv, dr, nn = map(np.array, (dv, dr, nn))
res = np.vstack([dv - X @ np.linalg.lstsq(X, dv, rcond=None)[0], dr - X @ np.linalg.lstsq(X, dr, rcond=None)[0]]).T
dof = len(dv) - X.shape[1]; cov = res.T @ res / dof
out['noise'] = dict(source='R120 (used): residuals of data - law after OLS on [1, I(2..13)]; sd per L-function at the mean number of spacings',
                    n_ref=float(nn.mean()), sd_var=math.sqrt(cov[0, 0]), sd_r=math.sqrt(cov[1, 1]), corr=cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1]))
print('noise (R120, %.0f spacings): sd Var %.5f  sd r %.5f  corr %+.3f' % (nn.mean(), out['noise']['sd_var'], out['noise']['sd_r'], out['noise']['corr']))
s = json.dumps(out, indent=1, default=float)
open(H + '/predictions_spacing.json.new', 'w').write(s)
hsh = hashlib.sha256(s.encode()).hexdigest(); open(H + '/predictions_spacing.sha256', 'w').write(hsh + '  predictions_spacing.json\n')
print('sha256', hsh)
