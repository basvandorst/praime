"""consec_fulltable.py: chi^2 of the full label tables (p_n, p_(n+k)) against N1, N2-LOS and N2-emp.
  N1      independent Chebotarev labels;                                   dof A^2 - 1
  N2-LOS  LOS for the abelian sign (consec_null.py) x Chebotarev inside;   dof A^2 - 1 (fully specified)
  N2-emp  the measured abelian table x Chebotarev inside each fibre;       dof A^2 - (number of abelian cells)
usage: consec_fulltable.py DIR X0 X1"""
import sys, os, glob, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import chi2 as CHI2
from consec_null import tables, MARG, CHI
d, X0, X1 = sys.argv[1], int(float(sys.argv[2])), int(float(sys.argv[3]))
segs = [np.load(f) for f in sorted(glob.glob(os.path.join(d, 'seg_*.npz'))) if X0 <= int(np.load(f)['L']) < X1]
nul = tables(X0, X1); out = {}
for nm in CHI:
    if nm + '__X1' not in segs[0].files: continue
    P = np.array(MARG[nm]); ch = np.array(CHI[nm][0]); A = len(P)
    for k in (1, 2, 3):
        X = sum(s[nm + '__X%d' % k].astype(float) for s in segs); N = X.sum()
        E1 = N * np.outer(P, P); E2 = N * np.array(nul[nm]['N2'][str(k)])
        # N2-emp: abelian table from the data
        sv = sorted(set(ch.tolist())); Ab = {(s, t): X[np.ix_(ch == s, ch == t)].sum() for s in sv for t in sv}
        cond = np.array([P[i] / P[ch == ch[i]].sum() for i in range(A)])
        E3 = np.array([[Ab[(ch[i], ch[j])] * cond[i] * cond[j] for j in range(A)] for i in range(A)])
        res = {}
        for tag, E, dof in (('N1', E1, A * A - 1), ('N2_LOS', E2, A * A - 1), ('N2_emp', E3, A * A - len(Ab))):
            c2 = float(((X - E) ** 2 / E).sum()); res[tag] = dict(chi2=c2, dof=dof, p=float(CHI2.sf(c2, dof)))
        out['%s_k%d' % (nm, k)] = res
        print('%-6s k=%d  N1 chi2 %10.1f/%d   N2-LOS %8.1f/%d   N2-emp %6.1f/%d (p %.3f)' % (nm, k, res['N1']['chi2'], A * A - 1, res['N2_LOS']['chi2'], A * A - 1, res['N2_emp']['chi2'], A * A - len(Ab), res['N2_emp']['p']))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fulltable_%s.json' % os.path.basename(d.rstrip('/'))), 'w'), indent=1)
