// Part 5: the page's transformations (prime.html, read-only) with the Mersenne exponents below 1000 as sources.
// The whole page script runs in a vm context whose DOM is a universal stub; then transform(p, code, k, N, S, C) is
// called with the page defaults k = 2, start S = 1, circle C = 1000 (N = 1000), destinations outside the circle skipped
// (the stub's select value is '' so the 'wrap' branch is off). Baseline: the same rules with all 168 primes < 1000.
const fs = require('fs'), vm = require('vm');
const html = fs.readFileSync(require('path').join(__dirname, '../../riemann/prime.html'), 'utf8');
const js = html.match(/<script>\n([\s\S]*)<\/script>/)[1];
const mk = () => new Proxy(function () {}, {
  get(t, k) {
    if (k === Symbol.toPrimitive) return () => 0;
    if (k === Symbol.iterator) return function* () {};
    if (k === 'value') return '';
    if (k === 'checked' || k === 'hidden' || k === 'disabled') return false;
    if (k === 'then') return undefined;
    return mk();
  },
  apply() { return mk(); }, construct() { return mk(); }, set() { return true; }
});
const VAL = { zeros: '575', experiment: 'bit-left', sourceMode: 'primes', overflow: 'skip', limit: '1000', limitNum: '1000',
  parameter: '2', start: '1', chainMax: '50' };
const elem = id => new Proxy(function () {}, { get(t, k) { if (k === 'value') return VAL[id] ?? ''; if (k === Symbol.toPrimitive) return () => 0;
  if (k === Symbol.iterator) return function* () {}; if (k === 'checked' || k === 'hidden' || k === 'disabled') return false; if (k === 'then') return undefined; return mk(); },
  apply() { return mk(); }, construct() { return mk(); }, set() { return true; } });
const doc = new Proxy(function () {}, { get(t, k) { if (k === 'getElementById') return elem; if (k === Symbol.toPrimitive) return () => 0; if (k === 'then') return undefined; return mk(); } });
const ctx = { document: doc, window: mk(), navigator: mk(), localStorage: mk(), location: mk(), history: mk(),
  Worker: mk(), Blob: mk(), URL: mk(), performance: { now: () => Date.now() }, setTimeout: () => 0, clearTimeout: () => 0,
  setInterval: () => 0, requestAnimationFrame: () => 0, console, Math, BigInt, Number, String, JSON, Map, Set, Array,
  Object, Float64Array, Uint8Array, Int32Array, Uint32Array, Float32Array, Int8Array, Uint16Array, Int16Array, Promise,
  Symbol, Error, RegExp, Date, isFinite, parseInt, parseFloat, matchMedia: mk(), getComputedStyle: mk(),
  addEventListener: () => 0, ResizeObserver: mk(), MutationObserver: mk(), IntersectionObserver: mk(), TextEncoder, TextDecoder };
vm.createContext(ctx);
let err = null;
try { vm.runInContext(js + '\n;globalThis.__x={transform,experiments};', ctx, { timeout: 60000 }); }
catch (e) { err = e; }
if (!ctx.__x) { try { vm.runInContext('globalThis.__x={transform,experiments};', ctx); } catch (e) { console.error('cannot reach transform:', err, e); process.exit(1); } }
if (err) console.error("page script stopped early (stub DOM):", String(err.stack).slice(0, 600));
const { transform, experiments } = ctx.__x;
const MEXP = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279], E = new Set(MEXP);
const SRC = MEXP.filter(x => x < 1000);
const isP = n => { if (n < 2) return false; for (let d = 2; d * d <= n; d++) if (n % d === 0) return false; return true; };
const PR = []; for (let q = 2; q < 1000; q++) if (isP(q)) PR.push(q);
const k = 2, S = 1, C = 1000, N = S + C - 1;
const rows = [];
let totE = 0, totEh = 0, totP = 0, totPh = 0;
for (const [i, ex] of experiments.entries()) {
  const code = ex[2];
  const run = srcs => { let d = 0, h = 0, hs = [];
    for (const p of srcs) { let y; try { y = transform(p, code, k, N, S, C); } catch (e) { y = null; }
      if (y === null || y === p) continue; d++; if (E.has(y)) { h++; hs.push(p + '->' + y); } }
    return { d, h, hs }; };
  const a = run(SRC), b = run(PR);
  totE += a.d; totEh += a.h; totP += b.d; totPh += b.h;
  rows.push({ nr: i + 1, name: ex[0], code, dE: a.d, hE: a.h, hits: a.hs.join(' '), dP: b.d, hP: b.h });
}
console.log(JSON.stringify({ rows, totE, totEh, totP, totPh, nSrc: SRC.length, nPr: PR.length }));
