const fs=require("fs");
const html=fs.readFileSync(process.argv[2],"utf8");
const grab=(a,b)=>{const i=html.indexOf(a),j=html.indexOf(b,i);if(i<0||j<0)throw new Error(a);return html.slice(i,j);};
const S={N:30,v:"inh"};
const src=grab("const CHB=","function chCorn(")+grab("function chQM(","function chBig(")+"function gcd(a,b){while(b){const t=a%b;a=b;b=t;}return a;}\n"+grab("const INS={N:0};","function inParts(")+grab("const INB=","function inLine(");
const api=new Function("S",src+";return{inGet,inRounds,inNode,INB};")(S);
const out={};for(const N of[2,3,30,250,5000]){S.N=N;const F=api.inGet(),rows=[];for(let k=2;k<=N;k++)if(F.dep[k]>=0){const H=F.inf[k];rows.push(H?[k,F.dep[k],H.F1,H.Q1,H.F2,H.Q2,H.w1||[],H.c||0,H.w2||[],H.x]:[k,0]);}out[N]={rows,cnt:F.cnt,cum:F.cum,D:F.D,best:F.best};}
if(process.argv[4]){out.box=[];for(const line of fs.readFileSync(process.argv[4],"utf8").split("\n")){const t=line.trim();if(!t)continue;const n=BigInt(t);let r=-1,ok=false;api.INB.part=false;try{r=api.inRounds(n);const nd=api.inNode(n);ok=nd.ok;}catch(e){r=-2;}out.box.push([t,r,ok?1:0,api.INB.part?1:0]);}}
fs.writeFileSync(process.argv[3],JSON.stringify(out));
