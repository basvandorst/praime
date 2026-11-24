import json, glob, numpy as np
for f in sorted(glob.glob(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/within/data/comp_*.json')):
    d = json.load(open(f)); lam = d['lam']
    print('lam %.4f' % lam, 'tau1_y13', {q: round(v,5) for q,v in d['tau1_triple_y13'].items()}, 'tau1_y11', {q: round(v,5) for q,v in d['tau1_triple_y11'].items()})
    for q in ('3','5','7'):
        c = d[q]; ch = np.array(c['chain_k1'])
        print('   q=%s V_k1 %.4f V_k2 %.4f Bq %+.4f  chain tau1 %+.5f  sum_{j>=2} chain %+.5f  tail_k1 %s  tail_y11 %s tail_w3 %s tail_k2 %s  MQ tau1 %+.5f sumMQ %+.5f' % (q, c['V_k1'], c['V_k2'], c['Bq'], ch[1], ch[2:].sum(), {k: round(v,4) for k,v in c['tail_k1'].items()}, {k: round(v,4) for k,v in c['tail_k1_y11'].items()}, {k: round(v,4) for k,v in c['tail_k1_w3'].items()}, {k: round(v,4) for k,v in c['tail_k2'].items()}, c['MQ_k1'][1], sum(c['MQ_k1'][1:])))
