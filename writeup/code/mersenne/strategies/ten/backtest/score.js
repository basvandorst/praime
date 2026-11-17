const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve(process.argv[2]); const runs = JSON.parse(fs.readFileSync(process.argv[3], 'utf8')); const outFile = process.argv[4];
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const port = 9400 + Math.floor(Math.random() * 400);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', '--hide-scrollbars', `--remote-debugging-port=${port}`, '--window-size=900,1200', 'about:blank'], { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
const SCORE = `(()=>{const N=S.N,pos=S.pos,cols=S.cols,rows=S.rows,dw=Math.max(1,S.dw||1),gc=Math.ceil(cols/dw),gr=Math.ceil(rows/dw);
 const G=new Int32Array(gc*gr),cx=new Int32Array(N+1),cy=new Int32Array(N+1);let placed=0,coll=0;
 for(let n=1;n<=N;n++){if(pos[2*n]<0){cx[n]=-1;continue;}const x=Math.floor(pos[2*n]/dw),y=Math.floor(pos[2*n+1]/dw);cx[n]=x;cy[n]=y;placed++;if(G[y*gc+x])coll++;G[y*gc+x]=n;}
 if(!window.__lab||window.__lab.N!==N){const y0=Math.round(Math.pow(N,0.4027)),sm=[];for(let q=2;q<=y0;q++)if(!S.sieve[q])sm.push(q);const sv=new Uint8Array(N+1).fill(1);sv[0]=0;sv[1]=0;for(const q of sm)for(let m=2*q;m<=N;m+=q)sv[m]=0;let c=1;for(const q of sm)c*=1-1/q;
  const nb=(N>>12)+1,bp=new Float64Array(nb),bs=new Float64Array(nb),pb=new Float64Array(nb);for(let n=y0+1;n<=N;n++){if(!S.sieve[n])bp[n>>12]++;if(sv[n])bs[n>>12]++;}for(let i=0;i<nb;i++)pb[i]=bs[i]?bp[i]/bs[i]:0;const mk=seed=>{let a=seed>>>0;const rnd=()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};const d=new Uint8Array(N+1);let k=0,real=0;for(let n=2;n<=N;n++){if(!S.sieve[n])real++;if(n<=y0){d[n]=S.sieve[n]?0:1;k+=d[n];continue;}if(sv[n]&&rnd()<pb[n>>12]){d[n]=1;k++;}}return {d,k,real};};
  const decs=[];for(let i=0;i<8;i++)decs.push(mk(1000+i));let rr=1;for(let i=0;i<3;i++)rr*=1;const f=decs.reduce((u,v)=>u+v.k,0)/8/decs[0].real;window.__lab={N,y0,decs,ratio:f};}
 const L=window.__lab,P=new Uint8Array(N+1);for(let n=2;n<=N;n++)P[n]=S.sieve[n]?0:1;
 const offs=[];for(let dy=0;dy<=3;dy++)for(let dx=-3;dx<=3;dx++){if(dy===0&&dx<=0)continue;offs.push([dx,dy]);}
 const cnt=lab=>{const r=[];for(const [dx,dy] of offs){let np=0,a=0,b=0,pp=0;for(let n=1;n<=N;n++){if(cx[n]<0)continue;const x=cx[n]+dx,y=cy[n]+dy;if(x<0||x>=gc||y>=gr)continue;const m=G[y*gc+x];if(!m||G[cy[n]*gc+cx[n]]!==n)continue;np++;const u=lab[n],v=lab[m];a+=u;b+=v;pp+=u&v;}r.push([np,a,b,pp]);}return r;};
 const mi=([np,a,b,pp])=>{if(np<50)return 0;const t=[[pp,a,b],[a-pp,a,np-b],[b-pp,np-a,b],[np-a-b+pp,np-a,np-b]];let s=0;for(const [k,r,c] of t)if(k>0)s+=k/np*Math.log2(k*np/(r*c));return s;};
 const rp=cnt(P),rd=L.decs.map(d=>cnt(d.d));let st=0,stS=0,zz=0,nz=0,maxz=0,mo=null,pairs=0;
 for(let i=0;i<offs.length;i++){const m=mi(rp[i]);st+=m;pairs+=rp[i][0];let ms=0,mean=0,vs=0;for(const d of rd){ms+=mi(d[i]);mean+=d[i][3];}ms/=rd.length;mean/=rd.length;for(const d of rd)vs+=(d[i][3]-mean)**2;vs=vs/(rd.length-1);stS+=ms;if(rp[i][0]>=2000&&mean>=20&&rp[i][3]>=20){const R=q=>q[3]*q[0]/(q[1]*q[2]),Rp=R(rp[i]);let Rd=0;for(const d of rd)Rd+=R(d[i]);Rd/=rd.length;const z=(Rp-Rd)/(Rd*Math.sqrt(1/rp[i][3]+1/(rd.length*mean)));zz+=z*z;nz++;if(Math.abs(z)>Math.abs(maxz)){maxz=z;mo=offs[i];}}}
 const ex=S.exps.filter(e=>e<=N&&cx[e]>=0),pl=[];for(let n=2;n<=N;n++)if(P[n]&&cx[n]>=0)pl.push(n);
 const nn=set=>{let s=0;for(let i=0;i<set.length;i++){let b=1e18;for(let j=0;j<set.length;j++){if(i===j)continue;const d=(cx[set[i]]-cx[set[j]])**2+(cy[set[i]]-cy[set[j]])**2;if(d<b)b=d;}s+=Math.sqrt(b);}return s/set.length;};
 let ez=null;if(ex.length>=8&&pl.length>200){const ix=new Map(pl.map((p,i)=>[p,i])),real=nn(ex);let a=77;const rnd=()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};const ds=[];for(let t=0;t<200;t++){const set=ex.map(e=>{const i=ix.get(e),lo=Math.max(0,Math.floor(i*0.7)-3),hi=Math.min(pl.length-1,Math.ceil(i*1.3)+3);return pl[lo+Math.floor(rnd()*(hi-lo+1))];});ds.push(nn(set));}const mu=ds.reduce((u,v)=>u+v,0)/ds.length,sd=Math.sqrt(ds.reduce((u,v)=>u+(v-mu)**2,0)/(ds.length-1));ez=(real-mu)/sd;}
 return JSON.stringify({placed,coll,pairs,struct:1000*st,structSieve:1000*stS,beyond:nz?Math.sqrt(zz/nz):null,nz,maxz,mo,expz:ez,nex:ex.length,cols,rows,dw});})()`;
(async () => {
  let L; for (let i = 0; i < 50; i++) { try { L = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); break; } catch (e) { await sleep(200); } }
  const ws = new WebSocket(L.find(t => t.type === 'page').webSocketDebuggerUrl); let id = 0; const pend = new Map();
  ws.onmessage = m => { const d = JSON.parse(m.data); if (d.id && pend.has(d.id)) { pend.get(d.id)(d); pend.delete(d.id); } };
  await new Promise(r => ws.onopen = r);
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  const ev = async e => { const r = await send('Runtime.evaluate', { expression: e, returnByValue: true }); return r.result?.result?.value; };
  await send('Page.enable'); await send('Runtime.enable');
  await send('Page.navigate', { url: 'file://' + file + '?layout=rows&links=none' + (process.argv[5] ? '&n=' + process.argv[5] : '') }); await sleep(2500);
  const out = [];
  for (const r of runs) {
    const t0 = Date.now();
    await ev(`(()=>{S.layout=${JSON.stringify(r.layout)};${r.period ? `S.P=${r.period};` : 'S.P=12;'}${r.width ? `S.w=${r.width};` : 'S.w=sqW();'}S.grid=null;try{full();}catch(e){window.__err=String(e);}})()`);
    let ok = false; for (let i = 0; i < 300; i++) { await sleep(60); if (await ev(`S.layout===${JSON.stringify(r.layout)}&&!!S.grid`)) { ok = true; break; } }
    let sc = null; try { sc = JSON.parse(await ev(SCORE)); } catch (e) { sc = { err: String(e) }; }
    const tag = r.tag || r.layout + (r.period ? '_P' + r.period : '') + (r.width ? '_w' + r.width : '');
    out.push({ tag, ...r, ok, dt: (Date.now() - t0) / 1000, ...sc });
    console.log(tag, ok ? '' : 'FAILED', sc && sc.struct != null ? `struct ${sc.struct.toFixed(2)} sieve ${sc.structSieve.toFixed(2)} beyond ${sc.beyond == null ? '-' : sc.beyond.toFixed(2)} maxz ${sc.maxz.toFixed(1)} expz ${sc.expz == null ? '-' : sc.expz.toFixed(2)}` : JSON.stringify(sc));
  }
  fs.writeFileSync(outFile, JSON.stringify(out, null, 1)); ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
