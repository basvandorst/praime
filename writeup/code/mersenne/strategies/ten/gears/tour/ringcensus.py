import subprocess, sys, json
LIST=sys.argv[1]; RMAX=int(sys.argv[2])
res={}
for r in range(2,RMAX+1):
    out=subprocess.run(['./tour','list',LIST,'ring',str(r),'1'],capture_output=True,text=True).stdout.split()
    L=[int(x) for x in out]; bad=[n for n in L if (n*n-1)%r!=0]
    res[r]=L
    print("ring %2d: liars %5d, first %s, with n^2 != 1 mod r: %d %s"%(r,len(L),L[0] if L else '-',len(bad),bad[:4]),flush=True)
json.dump(res,open('ring_'+LIST.replace('.bin','')+'.json','w'))
best=[]
rs=sorted(res)
for i,a in enumerate(rs):
    for b in rs[i+1:]:
        c=set(res[a])&set(res[b])
        best.append((len(c),a+b,a,b,min(c) if c else None))
best.sort()
print("pairs of rings with the fewest common liars (count, total cells, r1, r2, first):",best[:12])
