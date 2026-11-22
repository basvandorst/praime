import json, math, numpy as np
R = json.load(open('/home/bas/praime/writeup/data/study/r39.json')); E = json.load(open('/home/bas/praime/writeup/data/study/r13e.json'))
rd = {w: R[w]['m']/R[w]['n'] for w in R if w in ('train','A','B','C')}
for w, ref in (('train','A'), ('A','train')):
    kap = (1 - R[w]['2']['P'][0] - R[w]['2']['per_row'][0]) / R[w]['2']['P'][0]
    rows = []
    for c in range(1, 21):
        if c < 20:
            d = R[w][str(2*c)]; Ps = d['P']; heur = np.mean(d['per_row']); Pr = R[ref][str(2*c)]['P']
        else:
            Ps = [p for p in R[w]['P_pooled40'] if p > 0]; Pr = [p for p in R[ref]['P_pooled40'] if p > 0]; heur = np.mean([1-p for p in Ps])
        null = np.mean([1-p for p in Ps])
        resp = np.mean([(1-P) - kap*P*(math.log(q/P)/math.log(rd[ref]/rd[w]))**2 for P, q in zip(Ps, Pr)])
        rows.append((2*c, E[w]['class_%d'%c][0], null, heur, resp))
    print(w, 'gap measured null heur(r39) linear-response')
    for r in rows: print('  %3d %.3f %.3f %.3f %.3f' % r)
    m = np.array(rows)
    for col, lab in ((2,'null'),(3,'heur'),(4,'resp')):
        print('  unweighted SSE over classes vs', lab, round(float(np.sum((m[:,1]-m[:,col])**2)),3))
