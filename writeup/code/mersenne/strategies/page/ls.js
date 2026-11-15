const LSN=["quadform","eisenstein","mersieve","gaptowers","polyspiral","digitrows","chebrace","totient","goldbach","dickman","primroot","pili","ord2","wieferich","cunningham"],LSW=["gaptowers"],LSP=["quadform","mersieve","polyspiral","digitrows","chebrace"];
function lsClassNo(D){let h=0,am=0;for(let a=1;3*a*a<=-D;a++)for(let b=-a+1;b<=a;b++){if((b*b-D)%(4*a))continue;const c=(b*b-D)/(4*a);if(c<a||(b<0&&a===c))continue;if(gcd(gcd(a,Math.abs(b)),c)===1){h++;if(b===0||b===a||a===c)am++;}}return [h,am];}
function lsLi(x){if(x<=2)return 0;const g=0.5772156649015329,L=Math.log(x);let s=g+Math.log(L),t=1;for(let k=1;k<200;k++){t*=L/k;const d=t/k;s+=d;if(d<1e-15*s)break;}return s-1.0451637801174928;}
function lsLiF(x){const L=Math.log(x);let s=0.5772156649015329+Math.log(L),t=1;for(let k=1;k<300;k++){t*=L/k;const d=t/k;s+=d;if(Math.abs(d)<1e-15*Math.abs(s)+1e-300)break;}return s;}
function lsRmLi(x){const mu=[0,1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0,1,1,-1,0];let s=0;for(let k=2;k<25;k++)if(mu[k]&&Math.pow(x,1/k)>1.0001)s+=mu[k]/k*lsLiF(Math.pow(x,1/k));return s;}
function lsPrimeTD(n){if(n<2)return false;if(n%2===0)return n===2;for(let d=3;d*d<=n;d+=2)if(n%d===0)return false;return true;}
function lsFFT(re,im,inv){const n=re.length;for(let i=1,j=0;i<n;i++){let b=n>>1;for(;j&b;b>>=1)j^=b;j^=b;if(i<j){let t=re[i];re[i]=re[j];re[j]=t;t=im[i];im[i]=im[j];im[j]=t;}}
 for(let len=2;len<=n;len<<=1){const a=2*Math.PI/len*(inv?-1:1),wr=Math.cos(a),wi=Math.sin(a);for(let i=0;i<n;i+=len){let cr=1,ci=0;for(let k=0;k<len/2;k++){const ur=re[i+k],ui=im[i+k],xr=re[i+k+len/2],xi=im[i+k+len/2],vr=xr*cr-xi*ci,vi=xr*ci+xi*cr;re[i+k]=ur+vr;im[i+k]=ui+vi;re[i+k+len/2]=ur-vr;im[i+k+len/2]=ui-vi;const nr=cr*wr-ci*wi;ci=cr*wi+ci*wr;cr=nr;}}}
 if(inv)for(let i=0;i<n;i++){re[i]/=n;im[i]/=n;}}
function LS(name,N,w,P,pos,comp){const out={cols:0,rows:0,dw:1,fx:null,fy:null};const isP=n=>n>1&&!comp[n];
 const hid=(Wd,Hd,f)=>{out.cols=Wd+4;out.rows=Hd+4;out.fx=new Float64Array(N+1);out.fy=new Float64Array(N+1);for(let n=1;n<=N;n++){const r=f(n);if(!r){out.fx[n]=NaN;out.fy[n]=NaN;pos[2*n]=-1;pos[2*n+1]=-1;continue;}out.fx[n]=2+r[0];out.fy[n]=2+r[1];pos[2*n]=Math.max(0,Math.min(out.cols-1,Math.floor(out.fx[n])));pos[2*n+1]=Math.max(0,Math.min(out.rows-1,Math.floor(out.fy[n])));}};
 const spf=()=>{const s=new Int32Array(N+1);for(let i=2;i<=N;i++)if(!s[i])for(let j=i;j<=N;j+=i)if(!s[j])s[j]=i;return s;};
 const walk=Y=>{let lo=Infinity,hi=-Infinity;for(let n=1;n<=N;n++){lo=Math.min(lo,Y[n]);hi=Math.max(hi,Y[n]);}const sp=Math.max(1,hi-lo);hid(1400,700,n=>[(n-1)/Math.max(1,N-1)*1400,(hi-Y[n])/sp*700]);out.walk={lo,hi};};
 if(name==="quadform"){const Q=Math.max(1,Math.min(1000,Math.round(P))),B=new Int32Array(N+1).fill(-1);for(let b=0;Q*b*b<=N;b++){const r=Q*b*b;for(let a=0;a*a+r<=N;a++){const n=a*a+r;if(n>=1&&B[n]<0)B[n]=b;}}
  out.dw=2;const k=1396/Math.sqrt(N),sq=Math.sqrt(Q),A=n=>Math.round(Math.sqrt(n-Q*B[n]*B[n]));hid(1400,1400,n=>B[n]<0?null:[A(n)*k,1396-B[n]*sq*k]);
  const M=4*Q,inC=new Int32Array(M),outC=new Int32Array(M);let drawn=0,tot=0;for(let p=2;p<=N;p++){if(!isP(p)||M%p===0)continue;tot++;if(B[p]>=0){drawn++;inC[p%M]++;}else outC[p%M]++;}
  let cin=0,mixed=0,ctot=0;for(let c=0;c<M;c++){if(inC[c]+outC[c]===0)continue;ctot++;if(inC[c]&&outC[c])mixed++;else if(inC[c])cin++;}
  const [h,am]=lsClassNo(-4*Q),def=0.5*(h+am-1)*lsLi(Math.sqrt(N))/lsLi(N);out.qf={Q,h,am,def,drawn,tot,cin,mixed,ctot};out.tip=n=>B[n]<0?` · not of the form a² + ${Q}·b²`:` · ${A(n)}² + ${Q===1?"":Q+"·"}${B[n]}²`;}
 else if(name==="eisenstein"){const B=new Int32Array(N+1).fill(-1);for(let b=0;3*b*b<=N;b++)for(let a=b;a*a+a*b+b*b<=N;a++){const n=a*a+a*b+b*b;if(n>=1&&B[n]<0)B[n]=b;}
  out.dw=2;const k=1396/Math.sqrt(N),A=n=>Math.round((Math.sqrt(4*n-3*B[n]*B[n])-B[n])/2);hid(1400,700,n=>B[n]<0?null:[(A(n)+B[n]/2)*k,698-B[n]*Math.sqrt(3)/2*k]);
  let d1=0,d2=0,t1=0,t2=0;for(let p=2;p<=N;p++){if(!isP(p))continue;if(p%3===1){t1++;if(B[p]>=0)d1++;}else if(p%3===2){t2++;if(B[p]>=0)d2++;}}
  out.ef={d1,d2,t1,t2,three:B[3]>=0};out.tip=n=>B[n]<0?" · not of the form a² + ab + b²":` · ${A(n)}² + ${A(n)}·${B[n]} + ${B[n]}²`;}
 else if(name==="mersieve"){const K=Math.max(1,Math.min(256,Math.round(P))),row=new Int32Array(N+1),cnt=new Int32Array(K+2);
  const mm=(a,b,q)=>q<67108864?(a*b)%q:((a*(b>>>16))%q*65536+a*(b&65535))%q,pw=(p,q)=>{let r=1,x=2%q,e=p;while(e>0){if(e&1)r=mm(r,x,q);x=mm(x,x,q);e=Math.floor(e/2);}return r;};
  let viol=0,eu=0,euOk=0;for(let p=2;p<=N;p++){if(!isP(p))continue;let r=K+1;for(let k=1;k<=K;k++){const q=2*k*p+1,m8=q%8;if(m8!==1&&m8!==7)continue;if(p<60&&q>=2**p-1)break;if(pw(p,q)===1){r=k;break;}}
   row[p]=r;cnt[r]++;if(r<=K&&(r%4===2||(r%4===1&&p%4!==3)||(r%4===3&&p%4!==1)))viol++;const e=p>3&&p%4===3&&lsPrimeTD(2*p+1);if(e)eu++;if(e===(r===1))euOk++;}
  let tot=0;for(let r=1;r<=K+1;r++)tot+=cnt[r];const st=1396/(K+1);hid(1400,1400,n=>isP(n)?[(n-1)/Math.max(1,N-1)*1396,(row[n]-1)*st+((n*0.6180339887498949)%1)*0.62*st]:null);
  out.ms={K,cnt:Array.from(cnt),tot,viol,eu,euOk};out.tip=n=>!isP(n)?"":row[n]<=K?` · 2^${n} − 1 has the factor ${(2*row[n]*n+1).toLocaleString("en-US")} (k = ${row[n]})`:` · no factor 2kp + 1 with k ≤ ${K}`;}
 else if(name==="gaptowers"){const pr=[];for(let n=2;n<=N;n++)if(isP(n))pr.push(n);const m=pr.length,W2=Math.max(2,Math.min(w,m)),bands=Math.ceil(m/W2),bh=new Int32Array(bands),base=new Int32Array(N+1),gc=new Int32Array(400);let gmax=0,gp=0;
  for(let i=0;i<m;i++){const g=(i+1<m?pr[i+1]:N+1)-pr[i],b=Math.floor(i/W2);if(g>bh[b])bh[b]=g;if(i+1<m){gc[Math.min(399,g)]++;if(g>gmax){gmax=g;gp=pr[i];}}}
  const top=new Int32Array(bands+1);for(let b=0;b<bands;b++)top[b+1]=top[b]+bh[b]+2;out.cols=W2;out.rows=top[bands];pos[2]=-1;pos[3]=-1;
  for(let i=0;i<m;i++){const b=Math.floor(i/W2),c=i%W2,end=i+1<m?pr[i+1]:N+1;for(let n=pr[i];n<end;n++){pos[2*n]=c;pos[2*n+1]=top[b]+bh[b]-1-(n-pr[i]);base[n]=pr[i];}}
  let ch=0;for(let g=1;g<400;g++)if(gc[g]>gc[ch])ch=g;let odd=0;for(let g=3;g<400;g+=2)odd+=gc[g];
  out.gt={m,W2,bands,ch,chc:gc[ch],twins:gc[2],gmax,gp,odd};out.tip=n=>n<2?"":` · gap after ${base[n].toLocaleString("en-US")}: ${(()=>{const i=pr.indexOf(base[n]);return (i+1<m?pr[i+1]:N+1)-base[n];})()}`;}
 else if(name==="polyspiral"){const Q=Math.max(3,Math.min(60,Math.round(P)));let M=0;while(1+Q*M*(M+1)/2<N)M++;const k=1396/Math.max(1,M),V=[];for(let i=0;i<=Q;i++){const t=Math.PI/2+2*Math.PI*i/Q;V.push([Math.cos(t),Math.sin(t)]);}
  const T=m=>1+Q*m*(m+1)/2,ring=n=>{let m=Math.max(1,Math.floor(Math.sqrt(2*(n-1)/Q))-1);while(T(m)<n)m++;while(m>1&&T(m-1)>=n)m--;return m;};
  out.dw=2;hid(2800,2800,n=>{if(n===1)return [1398,1398];const m=ring(n),r=n-T(m-1)-1,i=Math.floor(r/m),t=r%m+1,x=m*V[i][0]+t*(V[i+1][0]-V[i][0]),y=m*V[i][1]+t*(V[i+1][1]-V[i][1]);return [1398+k*x,1398-k*y];});
  let best=null;for(let j=1;j<=Q;j++){let o=0,e=0,c=0;for(let m=1;;m++){const f=Q*m*m/2+(j-Q/2)*m+1;if(f>N)break;c++;if(f>=2){e+=1/Math.log(f);if(isP(f))o++;}}if(!best||o/Math.max(e,1e-9)>best.r)best={j,o,e,c,r:o/Math.max(e,1e-9)};}
  out.ps={Q,M,best};out.tip=n=>{if(n===1)return " · centre";const m=ring(n),r=n-T(m-1)-1;return ` · ring ${m}, side ${Math.floor(r/m)+1}`;};}
 else if(name==="digitrows"){const Bs=Math.max(2,Math.min(36,Math.round(P))),ds=new Uint16Array(N+1);for(let n=1;n<=N;n++)ds[n]=ds[Math.floor(n/Bs)]+n%Bs;
  let K=0;const cnt=new Int32Array(2000),rk=new Int32Array(N+1),pc=new Int32Array(2000);for(let n=1;n<=N;n++){rk[n]=cnt[ds[n]]++;if(ds[n]>K)K=ds[n];if(isP(n))pc[ds[n]]++;}
  const f=(n,y)=>[(rk[n]+0.5)/cnt[ds[n]]*1396,y];hid(1400,1400,n=>f(n,(ds[n]-1)*1396/Math.max(1,K-1)));
  const empty=[],ruled=[],extra=[];for(let s=1;s<=K;s++){if(!cnt[s])continue;const sh=gcd(s,Bs-1)>1;if(sh)ruled.push(s);if(!pc[s]){empty.push(s);if(!sh)extra.push(s);}}
  const leakP=[];for(let n=2;n<=N&&leakP.length<12;n++)if(isP(n)&&gcd(ds[n],Bs-1)>1)leakP.push(n);out.dr={Bs,K,empty,ruled,extra:extra.filter(s=>s>1),leakP,bp:isP(Bs)};out.tip=n=>` · digit sum ${ds[n]} in base ${Bs}`;}
 else if(name==="chebrace"){const Q=Math.max(3,Math.min(1000,Math.round(P))),sqr=new Uint8Array(Q),unit=new Uint8Array(Q);for(let a=1;a<Q;a++)if(gcd(a,Q)===1){unit[a]=1;sqr[(a*a)%Q]=1;}
  let nR=0,nN=0,r1=0;for(let a=1;a<Q;a++)if(unit[a]){if(sqr[a])nR++;else nN++;if((a*a)%Q===1)r1++;}
  const Y=new Float64Array(N+1);let neg=0,cntp=0,first=0;for(let n=1;n<=N;n++){let d=0;if(isP(n)&&Q%n!==0&&unit[n%Q]){d=sqr[n%Q]?-nN:nR;cntp++;}Y[n]=Y[n-1]+d;if(d&&Y[n]<0){neg++;if(!first)first=n;}}
  walk(Y);out.series=Y;out.cr={Q,nR,nN,r1,neg,cntp,first,end:Y[N],pred:nN*lsLi(Math.sqrt(N))/2};out.tip=n=>` · ${n%Q} mod ${Q}`+(unit[n%Q]?(sqr[n%Q]?", a square":", not a square"):"");}
 else if(name==="goldbach"){let L=1;while(L<2*N+2)L<<=1;const re=new Float64Array(L),im=new Float64Array(L);for(let n=2;n<=N;n++)if(isP(n))re[n]=1;lsFFT(re,im,false);for(let i=0;i<L;i++){const a=re[i],b=im[i];re[i]=a*a-b*b;im[i]=2*a*b;}lsFFT(re,im,true);
  const g=new Int32Array(N+1);let gmax=1,gmin=1e9,gminN=0;for(let n=4;n<=N;n+=2){g[n]=Math.round((re[n]+(isP(n/2)?1:0))/2);if(g[n]>gmax)gmax=g[n];if(n>=1000&&g[n]<gmin){gmin=g[n];gminN=n;}}
  let zero=0;for(let n=4;n<=N;n+=2)if(g[n]===0)zero++;const band=[0,0,0,0,0,0,0,0];for(let n=6;n<=N;n+=2){const v=g[n]*Math.log(n)**2/n;if(n%3===0){band[0]+=v;band[1]++;}else{band[2]+=v;band[3]++;if(n%5===0){band[4]+=v;band[5]++;}else{band[6]+=v;band[7]++;}}}
  out.parity=new Uint8Array(N+1);for(let n=6;n<=N;n+=6)out.parity[n]=1;hid(1400,1400,n=>n%2||n<4?null:[(n-1)/Math.max(1,N-1)*1396,1396-g[n]/gmax*1396]);
  out.gb={gmax,gmin,gminN,zero,ratio:(band[0]/band[1])/(band[2]/band[3]),r5:(band[4]/band[5])/(band[6]/band[7]),ev:Math.floor(N/2)-1};out.tip=n=>n%2||n<4?"":` · ${g[n].toLocaleString("en-US")} ways as p + q`;}
 else if(name==="dickman"){const lp=new Int32Array(N+1);for(let p=2;p<=N;p++)if(isP(p))for(let j=p;j<=N;j+=p)lp[j]=p;
  let s2=0,s3=0;for(let n=2;n<=N;n++){if(lp[n]*lp[n]<=n)s2++;if(lp[n]**3<=n)s3++;}hid(1400,1400,n=>n<2?[0,1396]:[(n-1)/Math.max(1,N-1)*1396,(1-Math.log(lp[n])/Math.log(n))*1396]);
  out.dk={s2,s3,tot:N-1};out.tip=n=>n<2?"":` · largest prime factor ${lp[n].toLocaleString("en-US")}`;}
 else if(name==="primroot"){const s=spf(),gr=new Int32Array(N+1),mm=(a,b,q)=>q<67108864?(a*b)%q:((a*(b>>>16))%q*65536+a*(b&65535))%q,pw=(x,e,q)=>{let r=1;x%=q;while(e>0){if(e&1)r=mm(r,x,q);x=mm(x,x,q);e=Math.floor(e/2);}return r;};
  let G=1,Gp=0;const cnt=new Int32Array(1000);let mod8=0,sq=0;for(let p=3;p<=N;p++){if(!isP(p))continue;const fs=[];let m=p-1;while(m>1){const r=s[m];fs.push(r);while(m%r===0)m/=r;}
   let g=2;for(;g<p;g++){let ok=true;for(const r of fs)if(pw(g,(p-1)/r,p)===1){ok=false;break;}if(ok)break;}gr[p]=g;cnt[Math.min(999,g)]++;if(g>G){G=g;Gp=p;}if(g===2&&(p%8===1||p%8===7))mod8++;const rt=Math.round(Math.sqrt(g));if(rt*rt===g)sq++;}
  gr[2]=1;let tot=0;for(let g=1;g<1000;g++)tot+=cnt[g];const st=1396/G;hid(1400,1400,n=>!isP(n)?null:n===2?[0,0]:[(n-1)/Math.max(1,N-1)*1396,(gr[n]-1)*st+((n*0.6180339887498949)%1)*0.62*st]);
  const used=[];for(let g=2;g<=G;g++)if(cnt[g])used.push(g);out.pr={G,Gp,c2:cnt[2],c3:cnt[3],tot,mod8,sq,used:used.length};out.tip=n=>isP(n)&&n>2?` · least primitive root ${gr[n]}`:"";}
 else if(name==="pili"){const Y=new Float64Array(N+1);let pc=0,mx=-Infinity,mxn=0;for(let n=1;n<=N;n++){if(isP(n))pc++;Y[n]=n<2?0:pc-lsLiF(n);if(n>=2&&Y[n]>mx){mx=Y[n];mxn=n;}}
  let r=0,c=0;const a=Math.log(1000),b=Math.log(N);if(N>2000)for(let i=0;i<3000;i++){const x=Math.round(Math.exp(a+(b-a)*i/2999));r+=Y[x]/lsRmLi(x);c++;}
  Y[1]=Y[2];walk(Y);out.series=Y;out.pl={mx,mxn,end:Y[N],rl:lsRmLi(N),half:-0.5*lsLiF(Math.sqrt(N)),ratio:c?r/c:NaN};out.tip=n=>n<2?"":` · π(n) − li(n) = ${Y[n].toFixed(1)}`;}
 else if(name==="ord2"){const s=spf(),ind=new Int32Array(N+1),mm=(a,b,q)=>q<67108864?(a*b)%q:((a*(b>>>16))%q*65536+a*(b&65535))%q,pw=(x,e,q)=>{let r=1;x%=q;while(e>0){if(e&1)r=mm(r,x,q);x=mm(x,x,q);e=Math.floor(e/2);}return r;};
  out.parity=new Uint8Array(N+1);let mI=1,tot=0,pf=0,odd=0,i1=0,bad4=0;const ordOf=new Int32Array(N+1);for(let q=3;q<=N;q++){if(!isP(q))continue;let o=q-1,m=q-1;while(m>1){const r=s[m];while(m%r===0)m/=r;while(o%r===0&&pw(2,o/r,q)===1)o/=r;}
   ordOf[q]=o;ind[q]=(q-1)/o;tot++;if(ind[q]>mI)mI=ind[q];if(ind[q]===1)i1++;if(lsPrimeTD(o)){out.parity[q]=1;pf++;if(o>2&&(ind[q]%2||ind[q]%8===4))bad4++;}if(ind[q]%2===1&&(q%8===1||q%8===7))odd++;}
  out.dw=2;const L2=Math.log2(Math.max(2,mI));hid(1400,1400,n=>!isP(n)||n<3?null:[(n-1)/Math.max(1,N-1)*1396,Math.log2(ind[n])/L2*1396]);
  out.o2={tot,pf,odd,i1,mI,bad4};out.tip=n=>isP(n)&&n>2?` · 2 has order ${ordOf[n].toLocaleString("en-US")} mod ${n.toLocaleString("en-US")}, index ${ind[n].toLocaleString("en-US")}`+(out.parity[n]?`: a factor of 2^${ordOf[n]} − 1`:""):"";}
 else if(name==="wieferich"){const q=new Float64Array(N+1),E=new Set(MEXP);let wf=[],half=0,tot=0,best=1,bestp=0,ks=0;const ys=[];
  for(let p=3;p<=N;p++){if(!isP(p))continue;const P2=BigInt(p)*BigInt(p);let r=1n,x=2n,e=p-1;while(e>0){if(e&1)r=r*x%P2;x=x*x%P2;e=Math.floor(e/2);}const qq=Number((r-1n)/BigInt(p));q[p]=qq/p;tot++;ys.push(q[p]);if(qq===0)wf.push(p);else if(q[p]<best){best=q[p];bestp=p;}if(q[p]<0.5)half++;}
  ys.sort((a,b)=>a-b);for(let i=0;i<ys.length;i++)ks=Math.max(ks,Math.abs(ys[i]-i/ys.length),Math.abs(ys[i]-(i+1)/ys.length));
  out.dw=2;hid(1400,1400,n=>!isP(n)||n<3?null:[(n-1)/Math.max(1,N-1)*1396,q[n]*1396]);out.wf={wf,half,tot,best,bestp,ks};out.tip=n=>isP(n)&&n>2?` · 2^(p−1) ≡ 1 + ${Math.round(q[n]*n).toLocaleString("en-US")}·p (mod p²)`:"";}
 else if(name==="cunningham"){const pt=n=>n<=1000002?n>1&&!comp[n]:lsPrimeTD(n),len=new Int32Array(N+1);let mx=1,mxp=0,sg=0,eu=0,pred=0;const cnt=new Int32Array(20),E=MEXP.filter(k=>k<=N);
  for(let p=2;p<=N;p++){if(!isP(p))continue;let L=1,v=p;while(pt(2*v+1)){v=2*v+1;L++;}len[p]=L;cnt[L]++;if(L>1)sg++;if(L>mx){mx=L;mxp=p;}}
  for(const k of E)if(k>3&&k%4===3&&len[k]>1)eu++;for(let n=3;n<=N;n++)pred+=1/(Math.log(n)*Math.log(2*n+1));pred*=1.3203236316937391;
  const st=1396/mx;hid(1400,1400,n=>!isP(n)?null:[(n-1)/Math.max(1,N-1)*1396,(len[n]-1)*st+((n*0.6180339887498949)%1)*0.62*st]);
  const sgE=E.filter(k=>k>3&&len[k]>1),mxs=[];for(let p=2;p<=N&&mxs.length<8;p++)if(isP(p)&&len[p]===mx)mxs.push(p);out.cu={mx,mxp,mxs,sg,pred,cnt:Array.from(cnt.slice(1,mx+1)),eu,sgE};out.tip=n=>isP(n)?` · chain of ${len[n]} from here: ${n.toLocaleString("en-US")}`+(len[n]>1?` → ${(2*n+1).toLocaleString("en-US")}`+(len[n]>2?" → …":""):""):"";}
 else if(name==="totient"){const s=spf(),ph=new Float64Array(N+1);ph[1]=1;for(let n=2;n<=N;n++){const p=s[n],q=n/p;ph[n]=q%p===0?ph[q]*p:ph[q]*(p-1);}
  hid(1400,1400,n=>[(n-1)/Math.max(1,N-1)*1396,(1-ph[n]/n)*1396]);out.tip=n=>` · φ(n)/n = ${(ph[n]/n).toFixed(4)}`;}
 return out;}
