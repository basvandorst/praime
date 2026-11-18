import subprocess, sys, itertools, json, time
LIST=sys.argv[1]; DEG=[int(x) for x in sys.argv[2].split(',')]; C=int(sys.argv[3]); out=sys.argv[4]
res=[]; t0=time.time()
for d in DEG:
    for cs in itertools.product(range(-C,C+1),repeat=d):
        if cs[-1]==0: continue
        r=subprocess.run(['./tour','list',LIST,'poly',str(d),'1']+[str(c) for c in cs],capture_output=True,text=True)
        L=[int(x) for x in r.stdout.split()]
        res.append((len(L),L[0] if L else None,d,cs))
    print('degree',d,'done',round(time.time()-t0),'s',file=sys.stderr,flush=True)
json.dump(res,open(out,'w'))
for d in DEG:
    rd=[x for x in res if x[2]==d]; z=[x for x in rd if x[0]==0]
    rd.sort(key=lambda x:(x[0],-(x[1] or 10**30)))
    print("degree %d: %d wheels, %d with no liar on the list; best: %s"%(d,len(rd),len(z),[("x^%d = %s"%(d,x[3]),x[0],x[1]) for x in rd[:6]]))
    hist={}
    for x in rd: hist[min(x[0],9)]=hist.get(min(x[0],9),0)+1
    print("   liars -> number of wheels:",sorted(hist.items()))
