const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve(process.argv[2]); const specs = JSON.parse(fs.readFileSync(process.argv[3], 'utf8')); const outFile = process.argv[4];
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const port = 9400 + Math.floor(Math.random() * 400);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', `--remote-debugging-port=${port}`, '--window-size=900,1200', 'about:blank'], { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
const LIB = `window.GS=(()=>{const N=S.N,TAU=2*Math.PI,y0=Math.round(Math.pow(N,0.4027)),sm=[];for(let q=2;q<=y0;q++)if(!S.sieve[q])sm.push(q);const sv=new Uint8Array(N+1).fill(1);sv[0]=0;sv[1]=0;for(const q of sm)for(let m=2*q;m<=N;m+=q)sv[m]=0;
 const nb=(N>>12)+1,bp=new Float64Array(nb),bs=new Float64Array(nb),pb=new Float64Array(nb);for(let n=y0+1;n<=N;n++){if(!S.sieve[n])bp[n>>12]++;if(sv[n])bs[n>>12]++;}for(let i=0;i<nb;i++)pb[i]=bs[i]?bp[i]/bs[i]:0;
 const mk=seed=>{let a=seed>>>0;const rnd=()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};const d=new Uint8Array(N+1);for(let n=2;n<=N;n++){if(n<=y0){d[n]=S.sieve[n]?0:1;continue;}if(sv[n]&&rnd()<pb[n>>12])d[n]=1;}return d;};
 const decs=[];for(let i=0;i<8;i++)decs.push(mk(5000+i));const P=new Uint8Array(N+1);for(let n=2;n<=N;n++)P[n]=S.sieve[n]?0:1;
 const offs=[];for(let dy=0;dy<=3;dy++)for(let dx=-3;dx<=3;dx++){if(dy===0&&dx<=0)continue;offs.push([dx,dy]);}
 const map=(sp,n)=>{if(sp.t==="rows")return [n%sp.w,Math.floor(n/sp.w)];if(sp.t==="spiral"){const th=TAU*sp.c*Math.pow(n,sp.a),r=0.6*Math.sqrt(n);return [Math.round(r*Math.cos(th)),Math.round(r*Math.sin(th))];}if(sp.t==="phyl"){const th=TAU*sp.g*n,r=0.6*Math.sqrt(n);return [Math.round(r*Math.cos(th)),Math.round(r*Math.sin(th))];}if(sp.t==="quad"){const r=Math.floor(Math.sqrt(n/sp.k));return [n-sp.k*r*r,r];}const s=sp.s,x=n%s,y=Math.floor(n/s)%s,z=Math.floor(n/(s*s));return [x+sp.a*z,y+sp.b*z];};
 const mi=([np,a,b,pp])=>{if(np<50)return 0;const t=[[pp,a,b],[a-pp,a,np-b],[b-pp,np-a,b],[np-a-b+pp,np-a,np-b]];let s=0;for(const [k,r,c] of t)if(k>0)s+=k/np*Math.log2(k*np/(r*c));return s;};
 const score=(sp,lo,hi)=>{const cx=new Int32Array(hi+1),cy=new Int32Array(hi+1);let x0=1e9,x1=-1e9,y0b=1e9,y1=-1e9;for(let n=1;n<=hi;n++){const p=map(sp,n);cx[n]=p[0];cy[n]=p[1];if(p[0]<x0)x0=p[0];if(p[0]>x1)x1=p[0];if(p[1]<y0b)y0b=p[1];if(p[1]>y1)y1=p[1];}
  const gc=x1-x0+1,gr=y1-y0b+1;if(gc*gr>3e7)return null;const G=new Int32Array(gc*gr);for(let n=1;n<=hi;n++){cx[n]-=x0;cy[n]-=y0b;G[cy[n]*gc+cx[n]]=n;}
  const cnt=lab=>{const r=[];for(const [dx,dy] of offs){let np=0,a=0,b=0,pp=0;for(let n=lo+1;n<=hi;n++){const x=cx[n]+dx,y=cy[n]+dy;if(x<0||x>=gc||y>=gr)continue;const m=G[y*gc+x];if(m<=lo||G[cy[n]*gc+cx[n]]!==n)continue;np++;const u=lab[n],v=lab[m];a+=u;b+=v;pp+=u&v;}r.push([np,a,b,pp]);}return r;};
  const rp=cnt(P),rd=decs.map(d=>cnt(d));let st=0,stS=0,zz=0,nz=0,maxz=0,pairs=0;for(let i=0;i<offs.length;i++){st+=mi(rp[i]);pairs+=rp[i][0];let ms=0,mean=0;for(const d of rd){ms+=mi(d[i]);mean+=d[i][3];}ms/=rd.length;mean/=rd.length;stS+=ms;if(rp[i][0]>=2000&&mean>=20&&rp[i][3]>=20){const R=q=>q[3]*q[0]/(q[1]*q[2]);let Rd=0;for(const d of rd)Rd+=R(d[i]);Rd/=rd.length;const z=(R(rp[i])-Rd)/(Rd*Math.sqrt(1/rp[i][3]+1/(rd.length*mean)));zz+=z*z;nz++;if(Math.abs(z)>Math.abs(maxz))maxz=z;}}
  return {pairs,struct:1000*st,decoy:1000*stS,beyond:nz?Math.sqrt(zz/nz):null,nz,maxz};};
 return {score};})();"ok"`;
(async () => {
  let L; for (let i = 0; i < 50; i++) { try { L = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); break; } catch (e) { await sleep(200); } }
  const ws = new WebSocket(L.find(t => t.type === 'page').webSocketDebuggerUrl); let id = 0; const pend = new Map();
  ws.onmessage = m => { const d = JSON.parse(m.data); if (d.id && pend.has(d.id)) { pend.get(d.id)(d); pend.delete(d.id); } };
  await new Promise(r => ws.onopen = r);
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  const ev = async e => { const r = await send('Runtime.evaluate', { expression: e, returnByValue: true }); return r.result?.result?.value; };
  await send('Page.enable'); await send('Runtime.enable');
  await send('Page.navigate', { url: 'file://' + file + '?layout=rows&links=none&n=1000000' }); await sleep(2500);
  console.log('lib', await ev(LIB), await ev('S.N'));
  const out = [];
  for (const sp of specs) {
    const r = JSON.parse(await ev(`JSON.stringify({a:GS.score(${JSON.stringify(sp)},0,500000),b:GS.score(${JSON.stringify(sp)},500000,1000000)})`));
    out.push({ sp, in: r.a, out: r.b });
  }
  fs.writeFileSync(outFile, JSON.stringify(out)); console.log('done', out.length); ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
