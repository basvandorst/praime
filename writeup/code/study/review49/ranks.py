import json, cypari2
pari=cypari2.Pari(); pari.allocatemem(200*10**6)
FS=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']
out={}
for c in FS:
    E=pari.ellinit(c['ainv']); r=int(pari.ellanalyticrank(E)[0]); out[c['id']]={'r':r,'w':c['w'],'N':c['N']}
json.dump(out,open('ranks.json','w'))
from collections import Counter
print(Counter((v['w'],v['r']) for v in out.values()))
