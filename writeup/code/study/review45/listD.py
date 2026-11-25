import json, sys
from sympy import factorint
d=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/triple/predictions_kappa3.json'))
D=d['Fn']['D']
print(len(D))
rows=[]
for x in D:
    f=factorint(abs(x))
    rows.append((x, dict(f)))
for x,f in rows:
    print(x, f, end='; ')
print()
# patterns
pr=[x for x,f in rows if len(f)==1 and list(f.values())[0]==1]
print('prime |D|:', pr[:10], len(pr))
print('4||D:', [x for x,f in rows if f.get(2,0)==2][:10])
print('8|D:', [x for x,f in rows if f.get(2,0)==3][:10])
print('2,3,5,7 all:', [x for x,f in rows if all(p in f for p in (2,3,5,7))])
print('2*3:', [x for x,f in rows if 2 in f and 3 in f][:10])
print('3*5:', [x for x,f in rows if 3 in f and 5 in f][:10])
print('3,5,7:', [x for x,f in rows if 3 in f and 5 in f and 7 in f][:10])
