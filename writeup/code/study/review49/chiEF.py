import json, numpy as np, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
ids=sorted(i for i in A.LAW if i in A.S2 and i not in A.EXCL)
S=A.meas(ids); Tf=np.array([A.P[i]['full'] for i in ids])
X=np.array([[1.0]+[float(A.FS[i]['g2']==g) for g in ('0','2','b')]+[float(A.FS[i]['g3']==g) for g in ('0','2','3','b')] for i in ids])
XtXi=np.linalg.pinv(X.T@X); n,k=X.shape
cS=XtXi@X.T@S; res=S-X@cS; Creg=np.cov(res.T,ddof=k)*(n-1)/(n-k)
def chi(TT,rows):
    cT=XtXi@X.T@TT; R=S-TT; cR=XtXi@X.T@R; rr=R-X@cR; Cn=rr.T@rr/(n-k)
    d=(cS-cT)[rows].ravel(); Aa=XtXi[np.ix_(rows,rows)]
    return d@np.linalg.solve(np.kron(Aa,Creg),d), d@np.linalg.solve(np.kron(Aa,Cn),d)
for lab,TT in [('frozen law',Tf),('GUE',np.tile(np.array(A.GUE),(n,1))),('deg1',np.array([A.P[i]['deg1'] for i in ids]))]+[('EF_%s'%x,A.arr(ids,'full','D'+x)+A.arr(ids,'full','G'+x)) for x in ('100','300','1000','H')]:
    a=chi(TT,[1,2,3]); b=chi(TT,[4,5,6,7])
    print('%-11s chi2_p2 registered-cov %8.1f noise-cov %8.1f | chi2_p3 registered-cov %8.1f noise-cov %8.1f'%(lab,a[0],a[1],b[0],b[1]))
import glob
E={}
for f in glob.glob(A.D+'e2all/*.json'):
    r=json.load(open(f)); E[r['id']]=r
T=A.arr(ids,'full','full')
for Xs in ('1000','10000','100000'):
    sX=np.array([E[i]['X'][Xs]['s2X'] for i in ids]); D=A.arr(ids,'full','D'+Xs); G=A.arr(ids,'full','G'+Xs)
    Tc=T-D+sX; Ec=G+sX
    a=chi(Tc,[1,2,3]); b=chi(Tc,[4,5,6,7]); c=chi(Ec,[1,2,3]); d=chi(Ec,[4,5,6,7])
    print('X=%-6s corrected law: chi2_p2 noise-cov %.1f, chi2_p3 %.1f | lines(exact, finite band)+GUE: chi2_p2 %.1f, chi2_p3 %.1f'%(Xs,a[1],b[1],c[1],d[1]))
# which cells: law - EF_1000 contrast coefficients vs data se
cT=XtXi@X.T@Tf; cE=XtXi@X.T@(A.arr(ids,'full','D1000')+A.arr(ids,'full','G1000'))
R=S-Tf; cR=XtXi@X.T@R; rr=R-X@cR; Cn=rr.T@rr/(n-k); se=np.sqrt(np.outer(np.diag(XtXi),np.diag(Cn)))
pk=[3,7,11,15,19]
for j,g in zip((1,2,3),('|a2|=0','|a2|=2','2|N')):
    print('%-7s law-EF %s   data-law %s   se %s'%(g,A.fmt((cT-cE)[j][pk]),A.fmt((cS-cT)[j][pk]),A.fmt(se[j][pk],'%.4f')))
