import json, math
r = json.load(open('/home/bas/praime/writeup/data/study/r53g.json'))['compare']
for g, v in r.items():
    print(g, v['data'], 'o2 %+.4f%% o3 %+.4f%% step %+.4f%% rem/step %.3f  o3/sd %.1f' % (100*v['rel_o2'], 100*v['rel_o3'], 100*v['step'], v['remainder_over_step'], v['rel_o3']/v['sd_rel']))
reg = [72,100,150,200,250,300,350,400,454,460,482,500,520,538,556,580,610,640,670,700]
print('registered 20 gaps with |rel_o3|/sd > 3:', sum(abs(r[str(g)]['rel_o3'])/r[str(g)]['sd_rel'] > 3 for g in reg))
print('max rel_o3 g<=100 %.4f%%, g<=300 %.4f%%, g<=700 %.4f%% (at g=%s)' % tuple([100*max(r[str(g)]['rel_o3'] for g in map(int, r) if g <= m) for m in (100, 300, 700)] + [max(r, key=lambda k: r[k]['rel_o3'])]))
