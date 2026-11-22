import json, numpy as np, math
from scipy.stats import spearmanr, pearsonr
from scipy.special import digamma
R=json.load(open('m1009_1_1008.json')); q=1009
gm=np.array([r['gmin'] for r in R]); a=np.array([r['a'] for r in R])
arch=math.log(q/math.pi)+digamma((1+a)/2)
def pred(S):
    lhs=S-np.median(S)+1/(0.25+np.median(gm)**2)
    return np.sqrt(np.clip(1/np.maximum(lhs,1e-9)-0.25,0,None))
for name,S in [('full',2*np.array([r['ReLL'] for r in R])+arch)]+[(f'trunc {X}',2*np.array([r['tr'][X] for r in R])+arch) for X in ('30','199','999')]:
    p=pred(S); print(name,'spearman',round(spearmanr(p,gm)[0],3),'pearson 1/(1/4+g^2) vs S',round(pearsonr(1/(0.25+gm**2),S)[0],3))
