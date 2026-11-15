const fs=require('fs');const t=fs.readFileSync(process.argv[2]||'../../prime2.html','utf8');const js=t.split('<script>')[1].split('</script>')[0];
const pre=js.slice(js.indexOf('const LSN'),js.indexOf('function sqW'));
const lib=`function gcd(a,b){while(b){[a,b]=[b,a%b];}return a;}function sieve(n){const s=new Uint8Array(n+1);s[0]=1;s[1]=1;for(let i=2;i*i<=n;i++)if(!s[i])for(let j=i*i;j<=n;j+=i)s[j]=1;return s;}`;
const r=new Function(lib+pre+`const N=500000,comp=sieve(1000002),pos=new Int32Array(2*(N+1));const t0=Date.now();const o=LS("mersieve",N,708,1000,pos,comp);return [o.ms.cnt[1001],o.ms.viol,o.ms.euOk,(Date.now()-t0)/1000];`)();
console.log('K=1000: survivors',r[0],'violations',r[1],'Euler agrees',r[2],'time',r[3],'s');
