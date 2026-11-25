import numpy as np, json, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analysis as AN
from analysis import *
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
fl = json.load(open(SCR + 'fastlaw_out.json'))
for band in ('lo', 'hi'):
    dd = {int(k): v for k, v in json.load(open(SCR + f'data_Fn_{band}.json')).items()}
    Dsel = [D for D in Ds if D != 1153 and str(D) in fl]
    a = arrays(Dsel)
    a['k3'] = np.array([dd[D]['k3'] for D in Dsel])
    rr = np.array([fl[str(D)][band] for D in Dsel])           # kfP, sigma2, kQ, kQ11, kQall
    a['full'] = rr[:, 0] + rr[:, 2] + 3*rr[:, 1] - 2*ls
    a['Q0'] = a['full'] + rr[:, 4] - rr[:, 2]
    # diagonly / twopt for sub-bands not recomputed: use the full-band differences (small L-dependence)
    i = [Ds.index(D) for D in Dsel]
    a['diagonly'] = a['full'] + (np.array([F['diagonly'][j] for j in i]) - np.array([F['full'][j] for j in i]))
    a['twopt'] = a['full'] + (np.array([F['twopt'][j] for j in i]) - np.array([F['full'][j] for j in i]))
    r = readings(a); se = boot(a, 600, seed=11)
    z = r['m']/se[5:]
    print(f"band {band}: n={len(Dsel)} zeros/set {np.mean([dd[D]['nz'] for D in Dsel]):.0f}  T2 p={r['pT']:.3f}  max|m|={np.max(np.abs(r['m'])):.5f} max|m/se|={np.max(np.abs(z)):.2f} at l={ls[np.argmax(np.abs(z))]}")
    print(f"   gA={r['gA']:.4f}+-{se[0]:.4f} gQ={r['gQ']:.4f}+-{se[1]:.4f} goff={r['goff']:.4f}+-{se[2]:.4f}  chi2 p=2,3,5 vs law", [round(r['chi2'][p]['full'], 1) for p in (2, 3, 5)])
    print('   m x1e4', np.round(r['m']*1e4, 1).tolist())
