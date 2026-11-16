import numpy as np, zlib, struct
def sieve(n):
    s=np.ones(n+1,dtype=bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return s
def png(path,rgb):
    rgb=np.ascontiguousarray(rgb.astype(np.uint8)); h,w,_=rgb.shape
    raw=np.concatenate([np.zeros((h,1),np.uint8),rgb.reshape(h,w*3)],axis=1).tobytes()
    def ch(t,d): return struct.pack(">I",len(d))+t+d+struct.pack(">I",zlib.crc32(t+d)&0xffffffff)
    open(path,"wb").write(b"\x89PNG\r\n\x1a\n"+ch(b"IHDR",struct.pack(">IIBBBBB",w,h,8,2,0,0,0))+ch(b"IDAT",zlib.compress(raw,9))+ch(b"IEND",b""))
def _ramp(stops,t):
    t=np.clip(t,0,1); xs=np.array([s[0] for s in stops]); out=np.zeros(t.shape+(3,))
    for c in range(3): out[...,c]=np.interp(t,xs,[s[1][c] for s in stops])
    return out
SEQ=[(0,(14,24,40)),(0.35,(38,70,120)),(0.7,(111,163,214)),(1,(240,248,255))]
DIV=[(0,(56,189,248)),(0.5,(14,24,40)),(1,(244,114,182))]
HOT=[(0,(14,24,40)),(0.3,(90,40,110)),(0.6,(244,114,182)),(0.85,(253,230,138)),(1,(255,255,255))]
def seq(a,lo=None,hi=None,stops=SEQ):
    a=np.asarray(a,float); lo=np.nanmin(a) if lo is None else lo; hi=np.nanmax(a) if hi is None else hi
    return _ramp(stops,(a-lo)/max(hi-lo,1e-300))
def div(a,m=None):
    a=np.asarray(a,float); m=np.nanmax(np.abs(a)) if m is None else m
    return _ramp(DIV,0.5+0.5*a/m)
def up(img,k): return np.repeat(np.repeat(img,k,axis=0),k,axis=1)
def lines(series,w=1200,h=420,ylo=None,yhi=None,bg=(14,24,40),grid=()):
    img=np.zeros((h,w,3))+np.array(bg)
    allv=np.concatenate([np.asarray(s[1],float) for s in series]); ylo=allv.min() if ylo is None else ylo; yhi=allv.max() if yhi is None else yhi
    for gy in grid:
        r=int((yhi-gy)/(yhi-ylo)*(h-1))
        if 0<=r<h: img[r,:]=(38,54,76)
    for xs,ys,col,th in series:
        xs=np.asarray(xs,float); ys=np.asarray(ys,float); x0,x1=series[0][0][0],series[0][0][-1]
        px=(xs-x0)/(x1-x0)*(w-1); py=(yhi-ys)/(yhi-ylo)*(h-1)
        for i in range(len(xs)-1):
            n=int(max(abs(px[i+1]-px[i]),abs(py[i+1]-py[i])))+2
            X=np.linspace(px[i],px[i+1],n); Y=np.linspace(py[i],py[i+1],n)
            for d in range(th):
                xi=np.clip(X.astype(int),0,w-1); yi=np.clip(Y.astype(int)+d,0,h-1); img[yi,xi]=col
    return img
