import numpy as np, json
T=np.load('twins.npy'); W=T[T<2*10**9]
out={}
for k in (1,2,3):
    for name,arr,starts in (('in',W,len(W)-k),('ext',T,int(np.sum(T<2e9)))):
        d=arr[k:k+starts]-arr[:starts]
        d=d[d<=1200]
        v,c=np.unique(d,return_counts=True)
        out[f'{k}_{name}']={int(a):int(b) for a,b in zip(v,c)}
        print(k,name,int(c.sum()), 'non-mult-6:', int(c[v%6!=0].sum()))
json.dump(out,open('counts.json','w'))
