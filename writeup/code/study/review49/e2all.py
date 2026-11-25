import sys, os, json, time, numpy as np
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import e2
e2.XS=[1000,10**4,10**5]
part,nparts,budget=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3])
t0=time.time()
for c in e2.FS[part::nparts]:
    f='e2all/%s.json'%c['id']
    if os.path.exists(f): continue
    if time.time()-t0>budget: break
    t1=time.time(); r=e2.run(c); json.dump(r,open(f,'w')); print(c['id'],round(time.time()-t1,1),flush=True)
print('left',sum(not os.path.exists('e2all/%s.json'%c['id']) for c in e2.FS[part::nparts]))
