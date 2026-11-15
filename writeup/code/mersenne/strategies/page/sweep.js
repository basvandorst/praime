const fs=require('fs');const t=fs.readFileSync(process.argv[2]||'../../prime2.html','utf8');const js=t.split('<script>')[1].split('</script>')[0];
const pre=js.slice(js.indexOf('const LSN'),js.indexOf('function sqW'));
const lib=`function gcd(a,b){while(b){[a,b]=[b,a%b];}return a;}function sieve(n){const s=new Uint8Array(n+1);s[0]=1;s[1]=1;for(let i=2;i*i<=n;i++)if(!s[i])for(let j=i*i;j<=n;j+=i)s[j]=1;return s;}`;
const run=new Function(lib+pre+`
const N=500000,comp=sieve(1000002),res={quadform:[],mersieve:[],chebrace:[],polyspiral:[]};
for(let P=1;P<=100;P++){const pos=new Int32Array(2*(N+1));const o=LS("quadform",N,708,P,pos,comp);res.quadform.push(o.qf);}
for(const K of [16,64,256]){const pos=new Int32Array(2*(N+1));const o=LS("mersieve",N,708,K,pos,comp);res.mersieve.push(o.ms);}
for(let P=3;P<=40;P++){const pos=new Int32Array(2*(N+1));const o=LS("chebrace",N,708,P,pos,comp);res.chebrace.push(o.cr);}
for(let P=3;P<=40;P++){const pos=new Int32Array(2*(N+1));const o=LS("polyspiral",N,708,P,pos,comp);res.polyspiral.push(o.ps);}
return res;`);
const r=run();fs.writeFileSync('sweep.json',JSON.stringify(r));
console.log('quadform P, h, drawn/tot, 1/(2h), mixed classes:');console.log(r.quadform.map(q=>`${q.Q}:${q.h}:${(q.drawn/q.tot).toFixed(4)}:${(1/(2*q.h)).toFixed(4)}:${q.mixed}`).join(' '));
console.log('mersieve:',r.mersieve.map(m=>`K=${m.K} survivors ${m.cnt[m.K+1]} of ${m.tot} viol ${m.viol} euler ${m.euOk} row1 ${m.cnt[1]}`).join(' | '));
console.log('chebrace P:nR:nN:end:pred:neg:first');console.log(r.chebrace.map(c=>`${c.Q}:${c.nR}:${c.nN}:${c.end}:${Math.round(c.pred)}:${c.neg}:${c.first}`).join(' '));
console.log('polyspiral best corner ray P:j:obs:exp:ratio');console.log(r.polyspiral.map(p=>`${p.Q}:${p.best.j}:${p.best.o}:${p.best.e.toFixed(1)}:${p.best.r.toFixed(2)}`).join(' '));
