// Test of best.js: nextMersenneExponent finds exactly the Mersenne exponents below 12,000.
// 1. the chain from 1 must be 2, 3, 5, ..., 11213 (the 23 exponents below 12,000);
// 2. no prime in (11213, 12000] passes (so the next exponent is above 12,000);
// 3. for 40 primes in [4200, 12000] that reach the FFT Lucas-Lehmer test, the FFT residue equals the BigInt residue
//    (independent exact arithmetic), and every trial-factoring factor divides 2^q - 1 (checked with BigInt).
const B = require('./best.js');
const KNOWN = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213];
const t0 = performance.now(); let ok = true;
const chain = []; let p = 1;
while (p < 11213) { p = B.nextMersenneExponent(p); chain.push(p); }
const c1 = JSON.stringify(chain) === JSON.stringify(KNOWN);
console.log(`1. chain from 1: ${chain.join(' ')}  ${c1 ? 'PASS' : 'FAIL'}  (${((performance.now() - t0) / 1000).toFixed(1)} s)`); ok &&= c1;
let extra = [];
for (let q = 11214; q <= 12000; q++) if (B.isPrimeSmall(q) && B.testExponent(q).prime) extra.push(q);
console.log(`2. primes in (11213, 12000] passing: ${extra.length ? extra.join(' ') : 'none'}  ${extra.length ? 'FAIL' : 'PASS'}`); ok &&= !extra.length;
let nRes = 0, nBad = 0, nTF = 0;
for (let q = 4201; q <= 12000 && nRes < 40; q += 97) {
  let r = q; while (!B.isPrimeSmall(r)) r++;
  const f = B.trialFactor(r, B.tfDepth(r));
  if (f) { nTF++; const F = BigInt(f); let x = 1n, b = 2n, e = r; while (e) { if (e & 1) x = x * b % F; b = b * b % F; e >>= 1; } if (x !== 1n || f <= 1) nBad++; continue; }
  nRes++;
  if (B.llResidueFFT(r).residue !== B.llResidueBigInt(r)) { nBad++; console.log('  residue differs at', r); }
}
console.log(`3. ${nRes} FFT residues vs BigInt and ${nTF} TF factors checked: ${nBad} bad  ${nBad ? 'FAIL' : 'PASS'}`); ok &&= !nBad;
console.log(ok ? `ALL PASS (${((performance.now() - t0) / 1000).toFixed(1)} s)` : 'FAILED');
if (typeof process !== 'undefined') process.exitCode = ok ? 0 : 1;
