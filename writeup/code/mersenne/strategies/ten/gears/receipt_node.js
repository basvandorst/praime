const fs=require("fs"),vm=require("vm"),crypto=require("crypto");
const ctx={console,BigInt,Math,Map,Set,Uint32Array,Array,Object,Number,performance};vm.createContext(ctx);vm.runInContext(fs.readFileSync(__dirname+"/receipt.js","utf8")+";this.rc={rcSha,rcG,rcRun,rcBuild,rcCheck,rcFull,rcK,rcLL};",ctx);
const rc=ctx.rc,hex=h=>Array.from(h).map(x=>x.toString(16).padStart(8,"0")).join("");
let shaOk=0;for(const s of["","abc","3|7|0|2,1|0,4|5,6","x".repeat(55),"y".repeat(56),"z".repeat(200)]){if(hex(rc.rcSha(s))===crypto.createHash("sha256").update(s).digest("hex"))shaOk++;}
const mode=process.argv[2];
if(mode==="make"){const out=[];for(const a of process.argv[3].split(",")){const p=Number(a),k=rc.rcK(p),F=rc.rcFull(p,k);out.push({p,k,prime:F.prime,ok:F.ok,turns:F.turns,y:[F.y[0].toString(),F.y[1].toString()],mus:F.mus.map(m=>[m[0].toString(),m[1].toString()]),tt:Math.round(F.tt),tb:Math.round(F.tb),tc:Math.round(F.tc*10)/10,ll:rc.rcLL(p)});}
 fs.writeFileSync(process.argv[4],JSON.stringify({sha:shaOk,rows:out}));}
else{const inp=JSON.parse(fs.readFileSync(process.argv[3],"utf8")),res=[];for(const r of inp){const G=rc.rcG(r.p),C=rc.rcCheck(G,r.y.map(BigInt),r.mus.map(m=>m.map(BigInt)));res.push({p:r.p,ok:C.ok,want:r.want});}
 fs.writeFileSync(process.argv[4],JSON.stringify(res));}
