import json,math,statistics as st,re
a={r['tag']:r for r in json.load(open('scores_5e5.json'))}; b={r['tag']:r for r in json.load(open('scores_1e6.json'))}
CON=set(re.search(r"CON=\{(.*?)\}",open('analyse.py').read()).group(1).replace("'","").split(','))
tags=[t for t in a if t in b and a[t].get('struct') is not None and b[t].get('struct') is not None and a[t]['layout'] not in CON and a[t]['pairs']>200000 and b[t]['pairs']>200000]
def rank(v):
    o=sorted(range(len(v)),key=lambda i:v[i]); r=[0]*len(v)
    for k,i in enumerate(o): r[i]=k
    return r
def corr(x,y):
    mx,my=st.mean(x),st.mean(y); return sum((p-mx)*(q-my) for p,q in zip(x,y))/math.sqrt(sum((p-mx)**2 for p in x)*sum((q-my)**2 for q in y))
sa=[a[t]['struct'] for t in tags]; sb=[b[t]['struct'] for t in tags]; ba=[a[t]['beyond'] or 0 for t in tags]; bb=[b[t]['beyond'] or 0 for t in tags]
ea=[(a[t]['expz'],b[t]['expz']) for t in tags if a[t]['expz'] is not None and b[t]['expz'] is not None]
print("layouts scored at both sizes: %d"%len(tags))
print("structure: rank correlation between N = 5e5 and 1e6: %.3f; median ratio 1e6 / 5e5: %.2f"%(corr(rank(sa),rank(sb)),st.median(y/x for x,y in zip(sa,sb) if x>5)))
print("beyond: rank correlation %.2f; median %.2f and %.2f; maximum %.2f (%s) and %.2f (%s)"%(corr(rank(ba),rank(bb)),st.median(ba),st.median(bb),max(ba),tags[ba.index(max(ba))],max(bb),tags[bb.index(max(bb))]))
print("exponent clustering z: correlation %.2f; mean %.2f and %.2f, sd %.2f and %.2f"%(corr([x for x,y in ea],[y for x,y in ea]),st.mean(x for x,y in ea),st.mean(y for x,y in ea),st.pstdev([x for x,y in ea]),st.pstdev([y for x,y in ea])))
hi=sorted(tags,key=lambda t:-(a[t]['beyond'] or 0)-(b[t]['beyond'] or 0))[:10]
print("highest 'beyond' at both sizes:"," ; ".join("%s %.2f/%.2f"%(t,a[t]['beyond'] or 0,b[t]['beyond'] or 0) for t in hi))
print("most structure (5e5 / 1e6):"," ; ".join("%s %.0f/%.0f"%(t,a[t]['struct'],b[t]['struct']) for t in sorted(tags,key=lambda t:-a[t]['struct'])[:14]))
ratio=[a[t]['struct']/a[t]['structSieve'] for t in tags if a[t]['structSieve']>5]+[b[t]['struct']/b[t]['structSieve'] for t in tags if b[t]['structSieve']>5]
print("structure of the primes / decoys: median %.3f, 5%%..95%% %.3f..%.3f"%(st.median(ratio),sorted(ratio)[len(ratio)//20],sorted(ratio)[-len(ratio)//20]))
