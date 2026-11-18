"""Collect the logs of the tournament into tournament.json and tournament.txt."""
import re, json, glob, os, struct
T={}
def parts(tag):
    out=[]
    for p in (0,1):
        e=open('big_%s_%d.err'%(tag,p)).read(); assert e.strip().endswith('done'),(tag,p)
        o=open('big_%s_%d.out'%(tag,p)).read().split('\n')
        out.append((e,[x for x in o if x.strip()]))
    return out
P=parts('misfit'); tested=0; liars=0; rings={}; last=0
for e,o in P:
    m=re.search(r'^tested (\d+) liars (\d+)',e,re.M); tested+=int(m.group(1)); liars+=int(m.group(2)); assert len(o)==int(m.group(2))
    last=max(last,int(re.search(r'last number read (\d+)',e).group(1)))
    for r,k in re.findall(r'ring (\d+) used (\d+) times',e): rings[int(r)]=rings.get(int(r),0)+int(k)
assert sum(rings.values())==tested
T['list_count']=tested; T['list_last']=last
T['misfit']={'tested':tested,'liars':liars,'rings':{str(k):rings[k] for k in sorted(rings)},'min_ring':min(rings),'max_ring':max(rings)}
T['census']={}; T['pair']=None
for tag in sorted(set(re.match(r'big_(census\d+)_',os.path.basename(f)).group(1) for f in glob.glob('big_census*_0.err')))+['cheap']:
    if tag=='cheap':
        if not (os.path.exists('big_cheap.err') and open('big_cheap.err').read().strip().endswith('done')): continue
        P=[(open('big_cheap.err').read(),[x for x in open('big_cheap.out').read().split('\n') if x.strip()])]
        T['pair']=None
    else:
        try: P=parts(tag)
        except (AssertionError,FileNotFoundError): continue
    for e,o in P:
        assert not o,(tag,o[:3])
        for r,a,b in re.findall(r'census ring (\d+): tested (\d+) liars (\d+)',e):
            c=T['census'].setdefault(r,{'tested':0,'liars':0}); c['tested']+=int(a); c['liars']+=int(b)
        m=re.search(r'pair \(.*\): tested (\d+) liars (\d+) strong liars (\d+)',e)
        if m and tag!='cheap' and os.path.exists('big_cheap.err'): m=None
        if m:
            if T['pair'] is None: T['pair']={'tested':0,'liars':0}
            T['pair']['tested']+=int(m.group(1)); T['pair']['liars']+=int(m.group(2))
        m=re.search(r'two-cell wheel with the doubling wheel: tested (\d+) liars (\d+) no wheel found (\d+)',e)
        if m:
            d=T.setdefault('duo_list',{'tested':0,'liars':0,'nostep':0}); d['tested']+=int(m.group(1)); d['liars']+=int(m.group(2)); d['nostep']+=int(m.group(3))
            for c,k in re.findall(r'two-cell wheel c = (\d+) used (\d+) times',e): d.setdefault('steps',{}); d['steps'][c]=d['steps'].get(c,0)+int(k)
        for r,a in re.findall(r'ring (\d+): (\d+) passed the half ring',e):
            T['census'].setdefault(r,{'tested':0,'liars':0})['half_ring_passes']=T['census'][r].get('half_ring_passes',0)+int(a)
if os.path.exists('big_ring2_5.err') and open('big_ring2_5.err').read().strip().endswith('done'):
    m=re.search(r'ring 5: tested (\d+) skipped (\d+) liars (\d+) last (\d+)',open('big_ring2_5.err').read())
    T['ring2_5']={'tested':int(m.group(1)),'skipped':int(m.group(2)),'liars':int(m.group(3))}
    assert not open('big_ring2_5.out').read().strip()
S={}
def single(fn):
    t=open(fn).read(); m=re.search(r'tested (\d+) liars (\d+)',t); return int(m.group(1)),int(m.group(2)),[int(x) for x in t.split('\n') if x.strip().isdigit()]
if os.path.exists('single.done'):
    a=single('single_doubling_1e8.txt'); S['doubling']={'first':single('single_doubling_first.txt')[2][0],'below1e8':a[1]}
    a=single('single_fib_1e8.txt'); S['fib']={'first':single('single_fib_first.txt')[2][0],'below1e8':a[1]}
    a=single('single_perrin_1e8.txt'); S['perrin']={'first':a[2][0],'below1e8':a[1],'list':a[2]}
    a=single('single_ring5_1e8.txt'); S['ring5']={'first':a[2][0],'below1e8':a[1]}
    T['odd_composites_below_1e8']=a[0]
    d=json.load(open('ring_psp2_to1e10.json')); S['ring5']['below1e10']=len(d['5'])
    T['fit_liars']=sum(len(v) for v in d.values()); T['fit_liars_misfits']=sum(1 for r,v in d.items() for n in v if (n*n-1)%int(r))
    T['ring5_classes_1e10']={str(c):sum(1 for n in d['5'] if n%5==c) for c in (1,2,3,4)}
if os.path.exists('single_frob_1e9.txt') and open('single_frob_1e9.txt').read().strip().endswith('done'):
    a=single('single_frob_1e9.txt'); S['duo']={'tested1e9':a[0],'below1e9':a[1]}
if os.path.exists('duo_row.json'): T['duo_row250']=json.load(open('duo_row.json'))['row250']
T['single']=S
pp=[5,7,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49]
done=set(int(r) for r in T['census']); R=4
for q in pp:
    if q in done: R=q
    else: break
T['all_r_upto']=[x for x in pp if x>R][0]-1 if R>4 else 4
if os.path.exists('pascal_check.json'): T['row250']=json.load(open('pascal_check.json'))['row250']
json.dump(T,open('tournament.json','w'),indent=1)
f=lambda n:format(n,',')
L=["Pascal's wheel on the smallest ring the number does not fit: %s candidates (all base-2 pseudoprimes below 2^64, last %d), %d came home"%(f(tested),last,liars)]
L.append("ring sizes used: "+", ".join("%d cells %s (%.2f %%)"%(k,f(rings[k]),100*rings[k]/tested) for k in sorted(rings)))
for r in sorted(T['census'],key=int):
    c=T['census'][r]; L.append("ring of %s cells, every candidate that does not fit it: %s tested, %d came home%s"%(r,f(c['tested']),c['liars'],(" (%s passed the half ring with a sign)"%f(c['half_ring_passes'])) if 'half_ring_passes' in c else ""))
if 'ring2_5' in T: L.append("second program (ring2.c), ring of 5 cells: %s tested (%s fit the ring or share a factor and were skipped), %d came home"%(f(T['ring2_5']['tested']),f(T['ring2_5']['skipped']),T['ring2_5']['liars']))
if T['pair']: L.append("doubling wheel and Fibonacci wheel, n = 2 or 3 (mod 5): %s tested, %d came home"%(f(T['pair']['tested']),T['pair']['liars']))
if 'duo' in S: L.append("two-cell wheel alone, every odd composite below 1e9: %s tested, %d came home"%(f(S['duo']['tested1e9']),S['duo']['below1e9']))
if 'duo_list' in T: L.append("two-cell wheel on all candidates below 2^64 (so: together with the doubling wheel): %s tested, %d came home, %d squares with no step; largest step used %s"%(f(T['duo_list']['tested']),T['duo_list']['liars'],T['duo_list']['nostep'],max(int(c) for c in T['duo_list'].get('steps',{'0':0}))))
for k,v in S.items():
    if k=='duo': continue
    L.append("single wheel %s: first mistake %s, mistakes among the odd composites below 1e8: %s%s"%(k,f(v['first']),f(v['below1e8']),(", below 1e10: %s"%f(v['below1e10'])) if 'below1e10' in v else ""))
if 'fit_liars' in T: L.append("rings of 2 to 40 cells below 1e10: %s (number, ring) mistakes, of which %d do not fit the ring; ring 5 by class of n mod 5: %s"%(f(T['fit_liars']),T['fit_liars_misfits'],T['ring5_classes_1e10']))
open('tournament.txt','w').write("\n".join(L)+"\n"); print("\n".join(L))
