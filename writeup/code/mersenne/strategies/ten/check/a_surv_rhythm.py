import numpy as np
N=500000; dead=np.zeros(N+1,bool); dead[:2]=True
for q in (2,3,5,7): dead[q::q]=True
S=np.nonzero(~dead)[0]; col=np.arange(len(S))%48; row=np.arange(len(S))//48
print("each column one remainder mod 210:",all(len(set((S[col==c]%210).tolist()))==1 for c in range(48)))
for q in (11,13,17):
    bad=0
    for c in range(48):
        r=row[(col==c)&(S%q==0)]; bad+=int((np.diff(r)!=q).sum())
    print("multiples of",q,": every",q,"th dot of each column; exceptions",bad)
