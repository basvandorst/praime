import sys,os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math, json
from stats import *
W=[(3*10**11,'pr_3e11.npy'),(6*10**11,'pr_6e11.npy'),(9*10**11,'pr_9e11.npy')]
L=10**8; H=10**6
prs=[(x0,np.load(f)) for x0,f in W]
res={}
for q in [3,4,10,210]:
    vals=[S1_blocks(pr,x0,L,H,q)[0] for x0,pr in prs]
    allv=np.concatenate(vals)
    res['S1_%d'%q]=dict(per=[v.mean() for v in vals], pooled=np.mean([v.mean() for v in vals]), err=allv.std(ddof=1)/math.sqrt(len(allv)))
    print('S1 q=%d'%q, ['%.4f'%v.mean() for v in vals], 'pooled %.4f err %.4f'%(res['S1_%d'%q]['pooled'],res['S1_%d'%q]['err']))
for q,k in [(3,1),(3,2),(10,1),(10,2)]:
    phi=len(red(q)); out={}
    parts={'chi':[], 'col':[], 'win':[], 'win2':[]}
    for x0,pr in prs:
        chi,col,win,win2,Ph,X=S2_blocks(pr,x0,L,H,q,k)
        nb=len(chi)
        parts['chi'].append(chi/(phi*(phi-1))); parts['col'].append(col/(phi-1)); parts['win'].append(win/((phi-1)**2)); parts['win2'].append(win2/((phi-1)**2))
    for key in parts:
        per=[v.sum()/(nb-1) for v in parts[key]]
        allv=np.concatenate(parts[key])
        err=allv.std(ddof=1)/math.sqrt(len(allv))*nb/(nb-1)
        print('S2 (%d,%d) %-5s'%(q,k,key), ['%.4f'%v for v in per], 'pooled %.4f err %.4f'%(np.mean(per),err))
