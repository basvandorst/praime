const fs=require("fs");
const html=fs.readFileSync(process.argv[2],"utf8");
const i=html.indexOf("const SFS={N:0}"),j=html.indexOf("function eBox(",i);if(i<0||j<0)throw new Error("block");
const S={N:30,v:"sf"};
const api=new Function("S","$","document",html.slice(i,j)+";return{sfGet,cnGet,cnCount,cnExact,cnAt,dTick,eLim};")(S,()=>({}),{});
const out={};
for(const N of[2,3,30,250,1000,5000]){S.N=N;S.v="sf";const F=api.sfGet();const lit=[];for(let k=2;k<=N;k++)if(F.lit[k])lit.push(k);
 const dec=[];for(let k=2;k<=N;k++)dec.push(api.dTick(k));const lims=[];for(let t=0;t<=F.T+1;t++)lims.push(api.eLim(t));
 S.v="sq";const sq=[];for(let k=2;k<=N;k++)sq.push(api.dTick(k));const sql=[];for(let t=0;t<=13;t++)sql.push(api.eLim(t));
 S.v="cnt";const C=api.cnGet();const pr=[];for(let k=2;k<=N;k++)if(C.pr[k])pr.push(k);const ex=[];for(const t of[0,1,2,3,5,7,10,29,30,31,100,241,999,4999,5000,6000])ex.push([t,C.V[api.cnExact(t)],C.C[api.cnExact(t)],C.mv[Math.min(t,N)]]);
 const mid=Math.min(N,7),CT=Array.from(api.cnAt(mid));
 out[N]={sf:{T:F.T,g:F.g,nx:F.nx,fr:F.fr,nw:F.nw,cum:F.cum,lit,dec,lims},sq:{dec:sq,lims:sql},cn:{n:C.n,x:C.x,pr,moves:C.mv[N],pi:C.C[0],V:C.V,C:Array.from(C.C),ex,mid,CT}};}
out.count=[];for(const x of[2,3,10,100,900,12345,1000000,99999999,1000000000,10000000000])out.count.push([x,api.cnCount(x)]);
fs.writeFileSync(process.argv[3],JSON.stringify(out));
