const fs=require("fs"),vm=require("vm");
const data=JSON.parse(fs.readFileSync(__dirname+"/tour/classpoly.json","utf8")).map(r=>[String(r[0]),r[2]]);
const src=fs.readFileSync(__dirname+"/chain.js","utf8").replace("__CHD__",JSON.stringify(data));
const ctx={console,BigInt,Math,Map,Uint8Array,Array,Object,Number,Error,JSON,performance};vm.createContext(ctx);vm.runInContext(src+";this.chTest=chTest;",ctx);
const mode=process.argv[2],out=[];
const ser=o=>JSON.stringify(o,(k,v)=>typeof v==="bigint"?v.toString():v);
if(mode==="range"){const top=Number(process.argv[3]),stop=BigInt(process.argv[4]||"30"),O=JSON.parse(process.argv[5]||"{}");for(let n=2;n<=top;n++){const r=ctx.chTest(BigInt(n),stop,O);out.push([String(n),r.ok,r.und?1:0,r.steps]);}}
else{const stop=BigInt(process.argv[4]||"10000000000"),O=JSON.parse(process.argv[5]||"{}");for(const line of fs.readFileSync(process.argv[3],"utf8").split("\n")){const s=line.trim();if(!s)continue;const t=performance.now(),r=ctx.chTest(BigInt(s),stop,O);out.push([s,r.ok,r.und?1:0,r.steps,Math.round(performance.now()-t),r.comp||"",r.nodes]);}}
fs.writeFileSync(process.argv[6]||"/dev/stdout",ser(out));
