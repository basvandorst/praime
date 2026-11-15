// Checks the Mersenne engine exactly as it runs on the page (../../riemann/prime.html): the functions are cut out of the
// page's script, the same text the Web Workers receive.
//   node page_check.js > page_check.json   then   python3 page_check.py < page_check.json
// 1. transformation 142: the chain from 1, and nextMersenneExponent(p) for every start p = 2 .. 11212, against the list
//    (Node has no Web Worker, so this is the page's synchronous path; it reaches 11213);
// 2. mxTest(q) for every prime q < 12000: exactly the listed exponents; every factor it reports divides 2^q - 1 (BigInt);
// 3. Lucas-Lehmer residues of the FFT (q >= 4200) for 18 primes, written out for page_check.py to compare with gmpy2;
// 4. mxIsZero (the zero test of the verifier) against the exact residue, on real runs and on constructed values
//    (0 and 2^q - 1 in plain and in balanced digits, and nonzero values).
const fs = require('fs'), path = require('path');
const html = fs.readFileSync(path.join(__dirname, '../../riemann/prime.html'), 'utf8');
const js = html.match(/<script>\n([\s\S]*)<\/script>/)[1];
const grab = name => { const i = js.indexOf('function ' + name + '('); if (i < 0) throw new Error(name); let d = 0, j = js.indexOf('{', i);
  for (let k = j; k < js.length; k++) { if (js[k] === '{') d++; if (js[k] === '}') { d--; if (!d) return js.slice(i, k + 1); } } };
const i0 = js.indexOf('function mxIsPrimeSmall('), i1 = js.indexOf('function mersTest(');
const engine = js.slice(i0, i1);
const m = new Function(grab('prime') + 'const mersCache=new Map(),mersNextCache=new Map();' + engine + grab('mersTest') + grab('mxIsZero') + grab('mxProve') +
  js.slice(js.indexOf('const MERS_SYNC='), js.indexOf('function mersWorkers(')) + grab('mersWorkers') + grab('nextMersenneExponent') +
  ';return {nextMersenneExponent,mxTest,mxSource,mxIsZero,mxProve,mxGetPlan,mxFftLengthFor,mxSqrStep,mxDigitsToBigInt,mxLlResidueFFT};')();
const LIST = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213], E = new Set(LIST);
const isP = n => { if (n < 2) return false; for (let d = 2; d * d <= n; d++) if (n % d === 0) return false; return true; };
const log = s => console.error(s);
let t = Date.now(), p = 1, chain = [];
while (true) { const y = m.nextMersenneExponent(p, Infinity); if (y === null) break; chain.push(y); p = y; }
log(`1. chain from 1: ${chain.join(' ')}; equals the list: ${chain.join() === LIST.join()} (${((Date.now() - t) / 1000).toFixed(1)} s)`);
let bad = 0; for (let q = 2; q < 11213; q++) { const y = m.nextMersenneExponent(q, Infinity); let e = q + 1; while (!E.has(e)) e++; if (y !== e) bad++; }
log(`   every start 2..11212: ${bad} wrong; start 11213 gives ${m.nextMersenneExponent(11213, Infinity)} (no workers in Node)`);
t = Date.now(); let found = [], fac = 0, facBad = 0;
for (let q = 2; q < 12000; q++) { if (!isP(q)) continue; const [ok, f] = m.mxTest(q); if (ok) found.push(q); if (f) { fac++; if ((2n ** BigInt(q)) % BigInt(f) !== 1n || BigInt(f) >= 2n ** BigInt(q) - 1n) facBad++; } }
log(`2. mxTest on every prime q < 12000: ${found.length} exponents, equal to the list: ${found.join() === LIST.join()}; ${fac} factors reported, ${facBad} wrong (${((Date.now() - t) / 1000).toFixed(1)} s)`);
const QS = [4201, 4253, 4423, 5003, 6007, 7001, 9689, 9941, 10007, 11213, 12011, 15013, 19937, 20011, 21701, 23209, 25013, 30011];
const rows = QS.map(q => { const r = m.mxLlResidueFFT(q), R = r.residue; return { q, n: r.n, maxErr: r.maxErr, lo64: (R & ((1n << 64n) - 1n)).toString(), m61: (R % ((1n << 61n) - 1n)).toString() }; });
log(`3. FFT residues of ${QS.length} primes written for page_check.py; largest rounding error ${Math.max(...rows.map(r => r.maxErr)).toFixed(4)}`);
let zc = 0, zb = 0;
for (const q of [4253, 4271, 9689, 9697, 11213, 11239, 19937, 19949, 23209, 23227]) { const N = m.mxFftLengthFor(q), P = m.mxGetPlan(q, N), x = new Float64Array(N); x[0] = 4;
  for (let i = 0; i < q - 2; i++) m.mxSqrStep(P, x, 2, 1); zc++; if (m.mxIsZero(P, x) !== (m.mxDigitsToBigInt(P, x) === 0n)) zb++; }
for (const q of [4253, 9941, 23209]) { const N = m.mxFftLengthFor(q), P = m.mxGetPlan(q, N), M = (1n << BigInt(q)) - 1n;
  const toD = v => { const x = new Float64Array(N); for (let i = 0; i < N; i++) x[i] = Number((v >> BigInt(P.c[i])) & ((1n << BigInt(P.c[i + 1] - P.c[i])) - 1n)); return x; };
  const bal = x => { const y = Float64Array.from(x); for (let i = 0; i < N - 1; i++) { const b = P.base[i]; if (y[i] >= b / 2) { y[i] -= b; y[i + 1] += 1; } } return y; };
  for (const v of [0n, M, 1n, M - 1n, 2n ** BigInt(q - 1), 12345678901234567890n % M]) for (const x of [toD(v), bal(toD(v))]) { zc++; if (m.mxIsZero(P, x) !== (v % M === 0n)) zb++; } }
log(`4. zero test of the verifier: ${zc} cases, ${zb} wrong`);
console.log(JSON.stringify(rows));
