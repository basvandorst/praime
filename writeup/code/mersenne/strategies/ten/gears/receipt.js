const RCK=(()=>{const k=[];let n=2;const pr=x=>{for(let d=2;d*d<=x;d++)if(x%d===0)return false;return true;};while(k.length<64){if(pr(n))k.push(n);n++;}return k;})();
const RCC=RCK.map(p=>Math.floor((Math.cbrt(p)%1)*4294967296)>>>0),RCH=RCK.slice(0,8).map(p=>Math.floor((Math.sqrt(p)%1)*4294967296)>>>0);
function rcSha(str){const m=[];for(let i=0;i<str.length;i++)m.push(str.charCodeAt(i)&255);const L=m.length;m.push(128);while(m.length%64!==56)m.push(0);const bl=L*8;m.push(0,0,0,Math.floor(bl/4294967296)&255,(bl>>>24)&255,(bl>>>16)&255,(bl>>>8)&255,bl&255);
 const H=RCH.slice(),w=new Uint32Array(64),rr=(x,n)=>(x>>>n)|(x<<(32-n));
 for(let o=0;o<m.length;o+=64){for(let i=0;i<16;i++)w[i]=(m[o+4*i]<<24|m[o+4*i+1]<<16|m[o+4*i+2]<<8|m[o+4*i+3])>>>0;for(let i=16;i<64;i++){const a=w[i-15],b=w[i-2];w[i]=(w[i-16]+(rr(a,7)^rr(a,18)^(a>>>3))+w[i-7]+(rr(b,17)^rr(b,19)^(b>>>10)))>>>0;}
  let[a,b,c,d,e,f,g,h]=H;for(let i=0;i<64;i++){const t1=(h+(rr(e,6)^rr(e,11)^rr(e,25))+((e&f)^(~e&g))+RCC[i]+w[i])>>>0,t2=((rr(a,2)^rr(a,13)^rr(a,22))+((a&b)^(a&c)^(b&c)))>>>0;h=g;g=f;f=e;e=(d+t1)>>>0;d=c;c=b;b=a;a=(t1+t2)>>>0;}
  H[0]=(H[0]+a)>>>0;H[1]=(H[1]+b)>>>0;H[2]=(H[2]+c)>>>0;H[3]=(H[3]+d)>>>0;H[4]=(H[4]+e)>>>0;H[5]=(H[5]+f)>>>0;H[6]=(H[6]+g)>>>0;H[7]=(H[7]+h)>>>0;}
 return H;}
function rcG(p){const N=(1n<<BigInt(p))-1n,P=BigInt(p),red=v=>{v=(v&N)+(v>>P);v=(v&N)+(v>>P);return v>=N?v-N:v;};
 return{p,N,T:p-2,x0:[2n%N,1n],red,sq:v=>[red(2n*v[0]*v[0]+N-1n),red(2n*v[0]*v[1])],mul:(u,v)=>{const m0=u[0]*v[0],m1=u[1]*v[1];return[red(m0+3n*m1),red((u[0]+u[1])*(v[0]+v[1])+2n*N*N-m0-m1)];},ser:v=>v[0].toString(16)+","+v[1].toString(16)};}
function rcPow(G,v,e){let r=[1n,0n];const b=e.toString(2);for(let i=0;i<b.length;i++){r=G.sq(r);if(b[i]==="1")r=G.mul(r,v);}return r;}
function rcHalves(T,k){const hs=[];for(let i=0;i<k;i++){if(T%2)T++;T/=2;hs.push(T);}return hs;}
function rcNeed(hs){const idx=new Set();let sums=[0];for(const h of hs){for(const s of sums)idx.add(s+h);sums=sums.concat(sums.map(s=>s+h));}return idx;}
function rcCh(G,x,y,mu,T){const h=rcSha(G.p+"|"+T+"|0|"+G.ser(x)+"|"+G.ser(y)+"|"+G.ser(mu));return(BigInt(h[0])<<32n)|BigInt(h[1]);}
function rcRun(G,k){const hs=rcHalves(G.T,k),idx=rcNeed(hs);let top=G.T;for(const i of idx)if(i>top)top=i;const saved=new Map();let v=G.x0,y=G.x0;for(let i=1;i<=top;i++){v=G.sq(v);if(idx.has(i))saved.set(i,v);if(i===G.T)y=v;}return{y,saved,top,hs};}
function rcBuild(G,k,y,R){const hs=R.hs,rs=[],memo=new Map(),f=(j,m)=>{if(j===0)return R.saved.get(m);const key=j+"|"+m;let v=memo.get(key);if(!v){v=G.mul(rcPow(G,f(j-1,m),rs[j-1]),f(j-1,m+hs[j-1]));memo.set(key,v);}return v;};
 let x=G.x0,T=G.T;const mus=[];for(let i=0;i<k;i++){if(T%2){y=G.sq(y);T++;}const mu=f(i,hs[i]),r=rcCh(G,x,y,mu,T);rs.push(r);mus.push(mu);x=G.mul(rcPow(G,x,r),mu);y=G.mul(rcPow(G,mu,r),y);T/=2;}return mus;}
function rcCheck(G,y,mus){let x=G.x0,T=G.T;for(const mu of mus){if(T%2){y=G.sq(y);T++;}const r=rcCh(G,x,y,mu,T);x=G.mul(rcPow(G,x,r),mu);y=G.mul(rcPow(G,mu,r),y);T/=2;}for(let i=0;i<T;i++)x=G.sq(x);return{ok:x[0]===y[0]&&x[1]===y[1],turns:T};}
function rcK(p){const T=p-2;let best=1,bc=Infinity;for(let k=1;k<=10;k++){const hs=rcHalves(T,k),chk=2*hs[k-1]+448*k,bld=224*(Math.pow(2,k)-k-1);if(bld<=0.35*2*T+1&&chk<bc){bc=chk;best=k;}}return best;}
function rcLL(p){if(p===2)return true;const N=(1n<<BigInt(p))-1n,P=BigInt(p);let s=4n;for(let i=0;i<p-2;i++){s=s*s-2n;s=(s&N)+(s>>P);s=(s&N)+(s>>P);if(s>=N)s-=N;}return s===0n;}
function rcFull(p,k){const G=rcG(p),t0=performance.now(),R=rcRun(G,k),t1=performance.now(),mus=rcBuild(G,k,R.y,R),t2=performance.now(),C=rcCheck(G,R.y,mus),t3=performance.now();return{p,k,prime:R.y[0]===0n,y:R.y,mus,ok:C.ok,turns:C.turns,T:G.T,top:R.top,saved:R.saved.size,tt:t1-t0,tb:t2-t1,tc:t3-t2};}
