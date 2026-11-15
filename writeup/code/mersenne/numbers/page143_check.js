const fs=require('fs');const html=fs.readFileSync(require('path').join(__dirname,'../../riemann/prime.html'),'utf8');const js=html.match(/<script>\n([\s\S]*)<\/script>/)[1];
const grab=name=>{const i=js.indexOf('function '+name+'(');let d=0,j=js.indexOf('{',i);for(let k=j;k<js.length;k++){if(js[k]==='{')d++;if(js[k]==='}'){d--;if(!d)return js.slice(i,k+1);}}};
const eng=js.slice(js.indexOf('function mxIsPrimeSmall('),js.indexOf('function mersTest('));
const m=new Function(grab('prime')+'const mersCache=new Map(),mersWant=new Set(),MERS_SYNC=3000;let mersPool=null;'+eng+grab('mersTest')+grab('mxIsZero')+grab('mersWorkers').replace('function mersWorkers(','function mersWorkersX(')+'function mersWorkers(){return [];}'+
 js.slice(js.indexOf('const A366_TOP='),js.indexOf('function nextMersenneExponent('))+';return {a366Test,a366Exact,a366Fft,a366A,nextA366};')();
const direct=k=>{const M=(1n<<BigInt(k))-1n,mod=M*BigInt(k)*BigInt(k);let e=(1n<<BigInt(k))-2n,b=BigInt(k-1)%mod,r=1n;while(e>0n){if(e&1n)r=r*b%mod;b=b*b%mod;e>>=1n;}return r===1n%mod;};
let t=Date.now(),bad=0,terms=[];for(let k=2;k<=2000;k++){const a=m.a366Test(k),b=direct(k);if(a)terms.push(k);if(a!==b){bad++;console.log('MISMATCH',k,a,b);}}
console.log(`1. a366Test = definition for every k in 2..2000: ${bad} mismatches; terms ${terms.join(' ')} (${((Date.now()-t)/1000).toFixed(0)} s)`);
const psp=[341,561,645,1105,1387,1729,1905,2047,2465,2701,2821,3277,4033,4369,4371,4681,5461,6601,7957,8321,8481,8911,10261,10585,11305];
console.log('2. base-2 pseudoprimes: (A) holds for all:',psp.every(m.a366A),'; congruence holds for:',psp.filter(m.a366Test));
t=Date.now();let fb=0,nf=0;for(const k of [4253,4423,9689,9941,11213,4201,4211,4217,4219,4229,5003,6007,7001,8009,9001,10007,11003,11939]){if(!m.a366A(k))continue;nf++;const a=m.a366Fft(k),b=m.a366Exact(k);if(a!==b){fb++;console.log('FFT MISMATCH',k,a,b);}}
console.log(`3. FFT path = exact BigInt on ${nf} k from 4201 to 11939: ${fb} mismatches (${((Date.now()-t)/1000).toFixed(0)} s)`);
t=Date.now();let p=1,ch=[];while(true){const y=m.nextA366(p,12000);if(y===null)break;ch.push(y);p=y;}
console.log(`4. chain from 1 (sync, up to 12000): ${ch.join(' ')} (${((Date.now()-t)/1000).toFixed(0)} s)`);
t=Date.now();console.log(`5. a366Test(65537) = ${m.a366Test(65537)}, a366Test(65539) = ${m.a366Test(65539)} (${((Date.now()-t)/1000).toFixed(0)} s)`);
