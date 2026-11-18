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
for rr in (5,7):
    fn='big_ring2_%d.err'%rr
    if os.path.exists(fn) and open(fn).read().strip().endswith('done'):
        m=re.search(r'ring %d: tested (\d+) skipped (\d+) liars (\d+) last (\d+)'%rr,open(fn).read())
        T['ring2_%d'%rr]={'tested':int(m.group(1)),'skipped':int(m.group(2)),'liars':int(m.group(3))}
        assert not open('big_ring2_%d.out'%rr).read().strip()
if all(os.path.exists(x) and 'exit 0' in open(x).read() for x in ('direct_pascal_a.txt','direct_pascal_b.txt')):
    T['direct_pascal']={'tested':0,'liars':0}
    for x in ('direct_pascal_a.txt','direct_pascal_b.txt'):
        m=re.search(r'tested (\d+) liars (\d+)',open(x).read()); T['direct_pascal']['tested']+=int(m.group(1)); T['direct_pascal']['liars']+=int(m.group(2))
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
if 'duo' in S:
    S['duo']['bound_exp']=9; S['duo']['tested']=S['duo']['tested1e9']; S['duo']['liars']=S['duo']['below1e9']
    if os.path.exists('single_frob_1e10.txt') and open('single_frob_1e10.txt').read().strip().endswith('done'):
        a=single('single_frob_1e10.txt'); S['duo']['bound_exp']=10; S['duo']['tested']+=a[0]; S['duo']['liars']+=a[1]
        ab=['single_frob_2e10a.txt','single_frob_2e10b.txt']
        if all(os.path.exists(x) and open(x).read().strip().endswith('done') for x in ab):
            for x in ab:
                a=single(x); S['duo']['tested']+=a[0]; S['duo']['liars']+=a[1]
            S['duo']['bound']=2*10**10
    S['duo'].setdefault('bound',10**S['duo']['bound_exp'])
    two=S['duo']['bound']==2*10**10
    S['duo']['bound_md']='2 · 10^10' if two else '10^%d'%S['duo']['bound_exp']
    S['duo']['bound_html']='2 × 10<sup>10</sup>' if two else '10<sup>%d</sup>'%S['duo']['bound_exp']
    S['duo']['bound_plain']='2e10' if two else '1e%d'%S['duo']['bound_exp']
if os.path.exists('duo_row.json'): T['duo_row250']=json.load(open('duo_row.json'))['row250']
T['single']=S
def done2(fn): return os.path.exists(fn) and open(fn).read().strip().endswith('done')
T['ring3']={}; T['folded']={}
for tag in 'abcdef':
    fs=['big3_%s_%d.err'%(tag,p) for p in (0,1)]
    if all(done2(x) for x in fs):
        nread=0
        for x in fs:
            e=open(x).read(); assert not open(x[:-3]+'out').read().strip(),x; nread+=int(re.search(r'numbers read (\d+)',e).group(1))
            for r,a,b,c in re.findall(r'ring (\d+): tested (\d+), passed the first step (\d+), came home (\d+)',e):
                d=T['ring3'].setdefault(r,{'tested':0,'first_step':0,'liars':0}); d['tested']+=int(a); d['first_step']+=int(b); d['liars']+=int(c)
        assert nread==T['list_count'],(tag,nread)
    fs=['big4_%s_%d.err'%(tag,p) for p in (0,1)]
    if all(done2(x) for x in fs):
        nread=0
        for x in fs:
            e=open(x).read(); assert not open(x[:-3]+'out').read().strip(),x; nread+=int(re.search(r'numbers read (\d+)',e).group(1))
            for m in re.finditer(r'ring (\d+): tested (\d+);(.*?) no fold: (\d+) tested, (\d+) passed the half ring; came home on the full ring: (\d+)',e):
                d=T['folded'].setdefault(m.group(1),{'tested':0,'folds':{},'silent':[],'nofold':0,'nofold_pass':0,'liars':0})
                d['tested']+=int(m.group(2)); d['nofold']+=int(m.group(4)); d['nofold_pass']+=int(m.group(5)); d['liars']+=int(m.group(6))
                for ee,a,b in re.findall(r'fold of (\d+) cells: (\d+) tested, (\d+) passed',m.group(3)):
                    q=d['folds'].setdefault(ee,{'tested':0,'passed':0}); q['tested']+=int(a); q['passed']+=int(b)
                d['silent']=sorted(set(d['silent'])|set(re.findall(r'fold of (\d+) cells: says nothing',m.group(3))),key=int)
        assert nread==T['list_count'],(tag,nread)
for r in T['folded']:
    if r in T['census']: assert T['folded'][r]['tested']==T['census'][r]['tested'],r
    if r in T['ring3']: assert T['folded'][r]['tested']==T['ring3'][r]['tested'],r
T['fold_pass_total']=sum(q['passed'] for d in T['folded'].values() for q in d['folds'].values())
pp=[5,7,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49]
done=set(int(r) for r in T['census'])|set(int(r) for r in T['ring3'])|set(int(r) for r in T['folded']); R=4
T['rings_done']=sorted(done)
assert all(d['liars']==0 for d in list(T['ring3'].values())+list(T['folded'].values()))
for q in pp:
    if q in done: R=q
    else: break
T['all_r_upto']=[x for x in pp if x>R][0]-1 if R>4 else 4
def _parts(m):
    out=[]; q=2
    while m>1:
        if m%q==0:
            k=1
            while m%q==0: m//=q; k*=q
            out.append((q,k))
        q+=1
    return out
def _settled(m):
    for q,k in _parts(m):
        lo={2:16,3:9}.get(q,q); x=lo
        while x<=k:
            if x not in done: return False
            x*=q
    return True
E=2
while _settled(E+2): E+=2
T['even_upto']=E; N=4
while _settled(N+1): N+=1
assert N==T['all_r_upto'],(N,T['all_r_upto'])
if os.path.exists('pascal_check.json'): T['row250']=json.load(open('pascal_check.json'))['row250']
json.dump(T,open('tournament.json','w'),indent=1)
f=lambda n:format(n,',')
L=["Pascal's wheel on the smallest ring the number does not fit: %s candidates (all base-2 pseudoprimes below 2^64, last %d), %d came home"%(f(tested),last,liars)]
L.append("ring sizes used: "+", ".join("%d cells %s (%.2f %%)"%(k,f(rings[k]),100*rings[k]/tested) for k in sorted(rings)))
for r in sorted(T['census'],key=int):
    c=T['census'][r]; L.append("ring of %s cells, every candidate that does not fit it: %s tested, %d came home%s"%(r,f(c['tested']),c['liars'],(" (%s passed the half ring with a sign)"%f(c['half_ring_passes'])) if 'half_ring_passes' in c else ""))
LOW={'25':5,'27':9,'32':16}
for r in sorted(T['ring3'],key=int):
    c=T['ring3'][r]; L.append("ring of %s cells (ring3.c), %s: %s tested, %s passed the first step (%s), %d came home"%(r,("every candidate that fits %d cells and not %s"%(LOW[r],r)) if r in LOW else "every candidate that does not fit it",f(c['tested']),f(c['first_step']),"the half ring with a sign" if r=='32' else "the half-size ring of the squares",c['liars']))
for r in sorted(T['folded'],key=int):
    c=T['folded'][r]; L.append("ring of %s cells (fold.c), every candidate that does not fit it: %s tested; %s%s%s; %d came home"%(r,f(c['tested']),"; ".join("fold onto %s cells: %s tested, %s passed"%(e,f(q['tested']),f(q['passed'])) for e,q in sorted(c['folds'].items(),key=lambda x:int(x[0]))),("; fold onto %s cells says nothing"%" and ".join(c['silent'])) if c['silent'] else "",("; no fold applies: %s tested on the half-size ring, %s passed"%(f(c['nofold']),f(c['nofold_pass']))) if c['nofold'] else "",c['liars']))
if T['folded']: L.append("candidates that passed a folded wheel, all rings together: %d"%T['fold_pass_total'])
L.append("every ring settled below 2^64 for numbers that do not fit it, without a gap, up to %d cells"%T['all_r_upto'])
for rr in (5,7):
    if 'ring2_%d'%rr in T: L.append("second program (ring2.c), ring of %d cells: %s tested (%s fit the ring or share a factor and were skipped), %d came home"%(rr,f(T['ring2_%d'%rr]['tested']),f(T['ring2_%d'%rr]['skipped']),T['ring2_%d'%rr]['liars']))
if 'direct_pascal' in T: L.append("Pascal's wheel on the smallest misfit ring, directly on every odd composite below 2e8 (no list, no shortcut): %s tested, %d came home"%(f(T['direct_pascal']['tested']),T['direct_pascal']['liars']))
if T['pair']: L.append("doubling wheel and Fibonacci wheel, n = 2 or 3 (mod 5): %s tested, %d came home"%(f(T['pair']['tested']),T['pair']['liars']))
if 'duo' in S: L.append("two-cell wheel alone, every odd composite below %s: %s tested, %d came home"%(S['duo']['bound_plain'],f(S['duo']['tested']),S['duo']['liars']))
if 'duo_list' in T: L.append("two-cell wheel on all candidates below 2^64 (so: together with the doubling wheel): %s tested, %d came home, %d squares with no step; largest step used %s"%(f(T['duo_list']['tested']),T['duo_list']['liars'],T['duo_list']['nostep'],max(int(c) for c in T['duo_list'].get('steps',{'0':0}))))
for k,v in S.items():
    if k=='duo': continue
    L.append("single wheel %s: first mistake %s, mistakes among the odd composites below 1e8: %s%s"%(k,f(v['first']),f(v['below1e8']),(", below 1e10: %s"%f(v['below1e10'])) if 'below1e10' in v else ""))
if 'fit_liars' in T: L.append("rings of 2 to 40 cells below 1e10: %s (number, ring) mistakes, of which %d do not fit the ring; ring 5 by class of n mod 5: %s"%(f(T['fit_liars']),T['fit_liars_misfits'],T['ring5_classes_1e10']))
open('tournament.txt','w').write("\n".join(L)+"\n"); print("\n".join(L))
