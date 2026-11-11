// Reads every prime up to N with transformation 139 of prime.html (exact window decoder,
// candidates coprime to 210, zeros up to 0.045 pi n), from the zeros in a file
// (for example the zeros recovered by primes_to_zeros.py), and checks each one.
// usage: node zeros_to_primes.js zeros.f64 N [depth]
const fs = require('fs'), path = require('path');
const html = fs.readFileSync(path.join(__dirname, 'prime.html'), 'utf8');
const js = html.match(/<script>\n([\s\S]*)<\/script>/)[1];
const a = js.indexOf('const Z30='), b = js.indexOf('let sieveCache');
const body = js.slice(a, b).replace('function zetaZeros(K){', 'function zetaZeros(K){if(globalThis.ZINJ)return globalThis.ZINJ;');
const buf = fs.readFileSync(process.argv[2]);
globalThis.ZINJ = Array.from(new Float64Array(buf.buffer, buf.byteOffset, buf.length / 8));
const m = new Function('$', body + ';return {riemannDecode};')(() => ({ value: 1576 }));
const K = ZINJ.length, N = Number(process.argv[3]), depth = Number(process.argv[4] || 4);
const isP = n => { if (n < 2) return false; if (n % 2 === 0) return n === 2; for (let q = 3; q * q <= n; q += 2) if (n % q === 0) return false; return true; };
const t = Date.now();
let p = 1, cnt = 0, bad = 0, last = 0;
for (;;) {
  const q = m.riemannDecode(p, 1, K, depth, true);
  if (q === null || q > N) break;
  let tq = p + 1; while (!isP(tq)) tq++;
  if (q !== tq) { bad++; if (bad < 4) console.log(`after ${p}: got ${q}, true ${tq}`); }
  cnt++; p = Math.max(q, tq); last = q;
}
let tot = 0; for (let n = 2; n <= N; n++) if (isP(n)) tot++;
console.log(`${K} zeros (the last at ${ZINJ[K - 1].toFixed(2)}): read ${cnt} primes up to ${last}, ${bad} wrong; primes up to ${N}: ${tot} [${((Date.now() - t) / 1000).toFixed(0)} s]`);
