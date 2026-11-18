const CHD=__CHD__.map(d=>[BigInt(d[0]),d[1].map(BigInt)]);
const CHB=100000,CHL=10000000000n;
const CHP=(()=>{const s=new Uint8Array(CHB+1),o=[];for(let i=2;i<=CHB;i++){if(!s[i]){o.push(i);for(let j=i*i;j<=CHB;j+=i)s[j]=1;}}return o;})();
const CHK=(()=>{const o=[];let ps=[],c=1;for(const p of CHP){if(c*p>=9007199254740992){o.push({c:BigInt(c),ps});ps=[];c=1;}ps.push(p);c*=p;}o.push({c:BigInt(c),ps});return o;})();
class ChC extends Error{}
function chPow(a,e,n){let r=1n%n;a%=n;while(e>0n){if(e&1n)r=r*a%n;a=a*a%n;e>>=1n;}return r;}
function chGcd(a,b){if(a<0n)a=-a;while(b){const t=a%b;a=b;b=t;}return a;}
function chSmall(n){if(n<2n)return false;const v=Number(n);for(const p of CHP){if(p*p>v)return true;if(v%p===0)return v===p;}return true;}
function chTrial(m){const fs=[];for(const ch of CHK){const p0=BigInt(ch.ps[0]);if(m<p0*p0)break;const r=Number(m%ch.c);for(const p of ch.ps)if(r%p===0){const P=BigInt(p);let e=0;while(m%P===0n){m/=P;e++;}fs.push([P,e]);}}
 if(m>1n&&m<CHL){fs.push([m,1]);m=1n;}return{fs,cof:m};}
function chSprp(n){if(n<3n||!(n&1n))return n===2n;let d=n-1n,s=0;while(!(d&1n)){d>>=1n;s++;}let x=chPow(2n,d,n);if(x===1n||x===n-1n)return true;for(let i=1;i<s;i++){x=x*x%n;if(x===n-1n)return true;}return false;}
function chJac(a,n){a%=n;if(a<0n)a+=n;let r=1;while(a){while(!(a&1n)){a>>=1n;const k=n&7n;if(k===3n||k===5n)r=-r;}const t=a;a=n;n=t;if((a&3n)===3n&&(n&3n)===3n)r=-r;a%=n;}return n===1n?r:0;}
function chIsqrt(n){if(n<2n)return n;let x=1n<<BigInt((n.toString(2).length+1)>>1);for(;;){const y=(x+n/x)>>1n;if(y>=x)return x;x=y;}}
function chInv(x,n){x%=n;if(x<0n)x+=n;let a=x,b=n,u=1n,v=0n;while(b){const q=a/b,t=a-q*b;a=b;b=t;const w=u-q*v;u=v;v=w;}if(a!==1n)throw new ChC("factor");u%=n;return u<0n?u+n:u;}
function chSqrtMod(a,p){a%=p;if(a<0n)a+=p;if(a===0n)return 0n;if(chPow(a,(p-1n)>>1n,p)!==1n)return null;let r;
 if((p&3n)===3n)r=chPow(a,(p+1n)>>2n,p);else{let q=p-1n,s=0;while(!(q&1n)){q>>=1n;s++;}let z=2n;while(chJac(z,p)!==-1){z++;if(z>2000n)return null;}let m=s,c=chPow(z,q,p),t=chPow(a,q,p);r=chPow(a,(q+1n)>>1n,p);
  while(t!==1n){let i=0,t2=t;while(t2!==1n){t2=t2*t2%p;i++;if(i===m)return null;}let b=c;for(let k=0;k<m-i-1;k++)b=b*b%p;m=i;c=b*b%p;t=t*c%p;r=r*b%p;}}
 return r*r%p===a?r:null;}
function chCorn(D,N){if(chJac(D,N)!==1)return null;let x0=chSqrtMod(D,N);if(x0===null)return null;if((x0-D)&1n)x0=N-x0;let a=2n*N,b=x0;const lim=chIsqrt(4n*N);while(b>lim){const t=a%b;a=b;b=t;}
 const r=4n*N-b*b,d=-D;if(r%d)return null;const c=r/d,s=chIsqrt(c);return s*s===c?[b,s]:null;}
function chTraces(D,t,s){if(D===-3n){const u=(t+3n*s)/2n,v=(t-3n*s)/2n;return[t,-t,u,-u,v,-v];}if(D===-4n)return[t,-t,2n*s,-2n*s];return[t,-t];}
function chAdd(P,Q,a,N){if(P===null)return Q;if(Q===null)return P;let num,den;if(P[0]===Q[0]){if((P[1]+Q[1])%N===0n)return null;if(P[1]!==Q[1])throw new ChC("roots");num=3n*P[0]*P[0]+a;den=2n*P[1];}else{num=Q[1]-P[1];den=Q[0]-P[0];}
 const l=num%N*chInv(den,N)%N;let x=(l*l-P[0]-Q[0])%N;if(x<0n)x+=N;let y=(l*(P[0]-x)-P[1])%N;if(y<0n)y+=N;return[x,y];}
function chMul(k,P,a,N){let R=null;const b=k.toString(2);for(let i=0;i<b.length;i++){R=chAdd(R,R,a,N);if(b[i]==="1")R=chAdd(R,P,a,N);}return R;}
function chPM(u,v,H,N){const h=H.length-1,r=new Array(2*h-1).fill(0n);for(let i=0;i<h;i++)if(u[i])for(let j=0;j<h;j++)r[i+j]=(r[i+j]+u[i]*v[j])%N;
 for(let i=2*h-2;i>=h;i--){const c=r[i];if(c)for(let j=0;j<=h;j++){const t=(r[i-h+j]-c*H[j])%N;r[i-h+j]=t<0n?t+N:t;}}return r.slice(0,h);}
function chPP(base,e,H,N){const h=H.length-1;let r=new Array(h).fill(0n);r[0]=1n;const b=e.toString(2);for(let i=0;i<b.length;i++){r=chPM(r,r,H,N);if(b[i]==="1")r=chPM(r,base,H,N);}return r;}
function chTrim(w,N){w=w.map(x=>{x%=N;return x<0n?x+N:x;});while(w.length&&w[w.length-1]===0n)w.pop();return w;}
function chPG(u,v,N){u=chTrim(u,N);v=chTrim(v,N);while(v.length){const iv=chInv(v[v.length-1],N);v=v.map(x=>x*iv%N);while(u.length>=v.length){const c=u[u.length-1],k=u.length-v.length;for(let j=0;j<v.length;j++)u[k+j]=u[k+j]-c*v[j];u=chTrim(u,N);}const t=u;u=v;v=t;}return u;}
function chRnd(R,N){let s="0x";const k=Math.ceil(N.toString(16).length/8)+1;for(let i=0;i<k;i++){R.s^=R.s<<13;R.s^=R.s>>>17;R.s^=R.s<<5;s+=(R.s>>>0).toString(16).padStart(8,"0");}return BigInt(s)%N;}
function chRoot(co,N,R){let H=co.slice().reverse().map(c=>{c%=N;return c<0n?c+N:c;}),h=H.length-1;if(h===1)return(N-H[0])%N;
 for(let it=0;it<40;it++){if(h===1)return(N-H[0])%N*chInv(H[1],N)%N;const base=new Array(h).fill(0n);base[0]=chRnd(R,N);base[1]=1n;const w=chPP(base,(N-1n)>>1n,H,N);w[0]=(w[0]-1n+N)%N;let g=chPG(H.slice(),w,N);
  if(g.length>1&&g.length<H.length){if(2*(g.length-1)>h){const q=[];let u=H.slice();while(u.length>=g.length){const c=u[u.length-1];q.push(c);const k=u.length-g.length;for(let j=0;j<g.length;j++){const t=(u[k+j]-c*g[j])%N;u[k+j]=t<0n?t+N:t;}u.pop();}g=q.reverse();const iv=chInv(g[g.length-1],N);g=g.map(x=>x*iv%N);}H=g;h=H.length-1;}}
 return null;}
function chCurve(D,co,m,q,N,R){let a0,b0;if(D!==-3n&&D!==-4n){const j=chRoot(co,N,R);if(j===null||j===0n||j===1728n%N)return null;const c=j*chInv(j-1728n,N)%N;a0=(N-3n*c%N)%N;b0=2n*c%N;}
 for(let it=0;it<60;it++){let x=chRnd(R,N),y,a,b;if(x===0n)continue;if(D===-3n){y=chRnd(R,N);a=0n;b=((y*y-x*x*x)%N+N)%N;}else if(D===-4n){y=chRnd(R,N);a=((y*y-x*x*x)%N+N)%N*chInv(x,N)%N;b=0n;}
  else{const z=(x*x*x+a0*x+b0)%N;if(z===0n)continue;a=a0*z%N*z%N;b=b0*z%N*z%N*z%N;x=z*x%N;y=z*z%N;}
  if(chGcd((4n*a*a*a+27n*b*b)%N,N)!==1n)continue;const Q=chMul(m/q,[x,y],a,N);if(Q===null)continue;if(chMul(q,Q,a,N)===null)return{t:"ec",N,D,a,b,x,y,m,q};}
 return null;}
function chQM(u,v,c,N){return[(u[0]*v[0]+c*u[1]*v[1])%N,(u[0]*v[1]+u[1]*v[0])%N];}
function chQP(u,e,c,N){let r=[1n,0n];const b=e.toString(2);for(let i=0;i<b.length;i++){r=chQM(r,r,c,N);if(b[i]==="1")r=chQM(r,u,c,N);}return r;}
function chMinus(N,fs){const out=[];for(const f of fs){let ok=false;for(let a=2n;a<300n;a++){if(chPow(a,N-1n,N)!==1n)throw new ChC("fermat");if(chGcd((chPow(a,(N-1n)/f,N)-1n+N)%N,N)===1n){out.push([f,a]);ok=true;break;}}if(!ok)return null;}return{t:"n-1",N,f:out};}
function chPlus(N,fs){let c=2n;for(;;){const j=chJac(c,N);if(j===-1)break;if(j===0){const g=chGcd(c,N);if(g!==1n&&g!==N)throw new ChC("factor");}c++;if(c>2000n)return null;}
 const out=[];for(const f of fs){let ok=false;for(let a=1n;a<300n;a++){const g=((a*a-c)%N+N)%N;if(chGcd(g,N)!==1n)continue;const gi=chInv(g,N),u2=chQM([a,1n],[a,1n],c,N),bt=[u2[0]*gi%N,u2[1]*gi%N];const e=chQP(bt,N+1n,c,N);if(e[0]!==1n||e[1]!==0n)throw new ChC("two cells");
   const w=chQP(bt,(N+1n)/f,c,N),nm=(((w[0]-1n)*(w[0]-1n)-c*w[1]*w[1])%N+N)%N;if(chGcd(nm,N)===1n){out.push([f,a]);ok=true;break;}}if(!ok)return null;}return{t:"n+1",N,c,f:out};}
function chRho(n,iters){let y=2n,q=1n,g=1n,r=1,done=0,x=y,ys=y;while(g===1n&&done<iters){x=y;for(let i=0;i<r;i++)y=(y*y+1n)%n;let k=0;while(k<r&&g===1n&&done<iters){ys=y;const lim=Math.min(100,r-k);for(let i=0;i<lim;i++){y=(y*y+1n)%n;q=q*(x>y?x-y:y-x)%n;}g=chGcd(q,n);k+=lim;done+=lim;}r*=2;}
 if(g===n){g=1n;for(let i=0;i<iters;i++){ys=(ys*ys+1n)%n;g=chGcd(x>ys?x-ys:ys-x,n);if(g>1n)break;}}return g>1n&&g<n?g:null;}
function chBig(c){for(let i=0;i<4;i++){if(chSprp(c))return c;const f=chRho(c,20000);if(f===null)return null;const o=c/f;c=f>o?f:o;}return chSprp(c)?c:null;}
function chCands(N,deep,O,menu){const cs=[],fits=q=>{const r=chIsqrt(q);return r>=2n&&(r-1n)**4n>N;};
 if(O.pm)for(const sg of[-1n,1n]){const kind=sg<0n?"n-1":"n+1",T=chTrial(N+sg);let F=1n;for(const[p,e]of T.fs)F*=p**BigInt(e);const big=G=>sg<0n?G*G>N:G>1n&&(G-1n)*(G-1n)>N;
   if(big(F)){if(menu)menu.push({kind,m:N+sg,q:1n});if(deep)continue;const part=[];let G=1n;for(const[p,e]of T.fs.slice().reverse()){part.push(p);G*=p**BigInt(e);if(big(G))break;}cs.push({q0:0n,kind,fs:part,q:null});}
   else if(T.cof>1n){let cof=T.cof;const pr=chSprp(cof);if(menu)menu.push({kind,m:N+sg,q:pr?cof:0n});if(deep===pr)continue;if(deep){cof=chBig(cof);if(cof===null||cof*cof<=4n*N)continue;}
    cs.push({q0:cof,kind,fs:big(cof)?[cof]:T.fs.map(f=>f[0]).concat([cof]),q:cof});}else if(menu)menu.push({kind,m:N+sg,q:0n});}
 if(O.quick&&cs.some(c=>c.q0===0n))return cs.filter(c=>c.q0===0n);
 if(O.ec&&N%2n&&N%3n)for(const[D,co]of CHD){if(co.length-1>O.h)continue;const ts=chCorn(D,N);if(ts===null)continue;for(const tr of chTraces(D,ts[0],ts[1])){const m=N+1n-tr,T=chTrial(m);let q=T.cof>1n?T.cof:(T.fs.length?T.fs[T.fs.length-1][0]:1n);
    const pr=T.cof<=1n||chSprp(T.cof),good=q>=2n&&fits(q);if(menu)menu.push({kind:"ec",D,m,q:good&&pr?q:0n});if(!good)continue;if(deep){if(pr)continue;q=chBig(T.cof);if(q===null||!fits(q))continue;}else if(!pr)continue;cs.push({q0:q,kind:"ec",D,co,m,q});}}
 const rank={"n-1":0,"n+1":1,ec:2};cs.sort((a,b)=>a.q0<b.q0?-1:a.q0>b.q0?1:rank[a.kind]-rank[b.kind]);return cs;}
function chProve(N,X){if(N<X.stop)return chSmall(N)?[{t:"td",N}]:null;const key=N.toString();if(X.memo.has(key))return X.memo.get(key);X.nodes++;let res=null;
 if(chSprp(N))for(const deep of[false,true]){if(res)break;if(deep&&!X.O.deep)break;const menu=!deep&&X.nodes===1?X.menu:null;for(const c of chCands(N,deep,X.O,menu)){if(X.nodes>X.max)break;let sub=[];if(c.q!==null&&c.q>=X.stop){sub=chProve(c.q,X);if(sub===null)continue;}
   const st=c.kind==="n-1"?chMinus(N,c.fs):c.kind==="n+1"?chPlus(N,c.fs):chCurve(c.D,c.co,c.m,c.q,N,X.R);if(st===null)continue;res=[st].concat(sub);break;}}
 X.memo.set(key,res);return res;}
function chTest(N,stop,O){O=Object.assign({pm:true,ec:true,h:4,deep:true},O||{});const X={stop:stop||CHL,memo:new Map(),nodes:0,max:4000,menu:[],O,R:{s:Number(N%1000003n)+7}};
 if(N<2n)return{ok:false,steps:null,menu:[],msg:"below 2"};if(N<X.stop)return{ok:chSmall(N),small:true,steps:chSmall(N)?[{t:"td",N}]:null,menu:[],nodes:0};if(!(N&1n))return{ok:false,even:true,steps:null,menu:[],nodes:0};
 if(!chSprp(N))return{ok:false,comp:"doubling",steps:null,menu:[],nodes:0};if(typeof n2Test==="function"&&!n2Test(N).ok)return{ok:false,comp:"norm2",steps:null,menu:[],nodes:0};
 try{const steps=chProve(N,X);return{ok:steps!==null,und:steps===null,steps,menu:X.menu,nodes:X.nodes};}catch(e){if(e instanceof ChC)return{ok:false,comp:e.message,steps:null,menu:X.menu,nodes:X.nodes};throw e;}}
