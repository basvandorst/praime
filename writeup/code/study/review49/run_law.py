import sys, json, time, os, numpy as np, math
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mylaw as M
FS=json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']
part,nparts,budget=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3])
BANDS={'full':(20,400,24),'lo':(20,200,16),'hi':(200,400,16),'b60':(60,400,20)}
XS=[10,30,100,300,1000,10**4,10**5,10**6,10**7]
t0=time.time()
for f in FS[part::nparts]:
    out='lawout/%s.json'%f['id']
    if os.path.exists(out): continue
    if time.time()-t0>budget: break
    t1=time.time(); cv=M.Curve(f['ainv'],f['N'])
    FJ=cv.FJ(); FK=cv.FK(); FKn=cv.FK('nosym'); FKd=cv.FK('deg1'); FJd=M.ZD+1/M.EPS**2
    good=~cv.bad
    def wdeg1(k):
        w=np.ones(len(M.PR)); w[cv.bad]=0.0; return w
    R={'id':f['id'],'N':f['N'],'LS1':cv.LS1}
    W2={k:cv.W(k)[0] for k in range(1,30)}; W3={k:cv.W(k)[1] for k in range(1,30)}
    bad2=bool(cv.bad[0]); bad3=bool(cv.bad[1])
    for b,(lo,hi,n) in BANDS.items():
        tn,wn=M.nodes(lo,hi,n); Ln=M.Lt(f['N'],tn); wt=wn*Ln; wt/=wt.sum()
        acc={}
        def add(k,v,w):
            acc[k]=acc.get(k,0)+w*np.asarray(v)
        for L,w,t in zip(Ln,wt,tn):
            corr=cv.diag_corr(L); base=M.GUE+corr
            add('full',base+M.s2_from_F(FJ,FK,L),w)
            add('nosym',base+M.s2_from_F(FJ,FKn,L),w)
            add('deg1',M.GUE+M.diag_corr_w(cv,L,wdeg1)+M.s2_from_F(FJd,FKd,L),w)
            z=np.zeros(20)
            if not bad2:
                add('dk1p2',sum(M.single(L,2,k,2-W2[k]) for k in range(2,30)),w)   # alt - full
                add('dw1p2',M.single(L,2,1,1-W2[1]),w)
                add('d2_2',M.single(L,2,2,W2[2]),w)
            else:
                add('dk1p2',z,w); add('dw1p2',z,w); add('d2_2',z,w)
            if not bad3:
                add('dk1p3',sum(M.single(L,3,k,2-W3[k]) for k in range(2,30)),w)
                add('dw1p3',M.single(L,3,1,1-W3[1]),w)
            else:
                add('dk1p3',z,w); add('dw1p3',z,w)
            dk1=np.zeros(20)
            for k in range(2,30):
                Wk=cv.W(k); sel=good&(M.PR.astype(float)**k<1e25)
                ph=np.pi*k*np.outer(M.LS,M.LNP[sel])/L
                dk1+=(2/np.pi**2)*(np.sin(ph)**2*((2-Wk[sel])/(k*k*M.PR[sel].astype(float)**k))).sum(1)
            add('dk1',dk1,w)
            for X in XS:
                add('D%d'%X,M.diag_X(cv,L,X),w); add('G%d'%X,M.gue_completion(L,X),w)
            XH=math.exp(L); add('DH',M.diag_X(cv,L,XH),w); add('GH',M.gue_completion(L,XH),w)
            XH2=math.exp(L/2); add('DH2',M.diag_X(cv,L,XH2),w); add('GH2',M.gue_completion(L,XH2),w)
            add('meanL',L,w)
        R[b]={k:(v.tolist() if hasattr(v,'tolist') else float(v)) for k,v in acc.items()}
    json.dump(R,open(out,'w'))
    print(f['id'], round(time.time()-t1,1), flush=True)
left=sum(not os.path.exists('lawout/%s.json'%f['id']) for f in FS[part::nparts])
print('left',left, round(time.time()-t0))
