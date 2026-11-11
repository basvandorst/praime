// Validates the zero formulas of prime.html against trial division.
// Each formula decodes circles that start at S with the true primes below S known
// (the page does the same), using just enough zeros for its height (times 1.02).
// usage: node validate.js [formulas]   formulas: 134,135,136,137,138,139,140 (default all)
const fs = require('fs'), path = require('path');
const html = fs.readFileSync(path.join(__dirname, 'prime.html'), 'utf8');
const js = html.match(/<script>\n([\s\S]*)<\/script>/)[1];
const a = js.indexOf('const Z30='), b = js.indexOf('let sieveCache');
const m = new Function('$', js.slice(a, b) + ';return {zetaZeros,riemannDecode,riemannClass,riemannMod24,riemannMod840,riemannAuto,decStates,clsStates,m24States,m840States,DEC_PRIMES};')(() => ({ value: 1576 }));
const isP = n => { if (n < 2) return false; if (n % 2 === 0) return n === 2; for (let q = 3; q * q <= n; q += 2) if (n % q === 0) return false; return true; };
const K = 60000, G = m.zetaZeros(K), P = Math.PI;
const cnt = T => { let lo = 0, hi = K; while (lo < hi) { const md = (lo + hi) >> 1; if (G[md] <= T) lo = md + 1; else hi = md; } return Math.min(K, lo + 1); };
const clear = () => { m.decStates.clear(); m.clsStates.clear(); m.m24States.clear(); m.m840States.clear(); };
const rho = k => m.DEC_PRIMES.slice(0, k).reduce((r, q) => r * (1 - 1 / q), 1);
const F = {
  134: { f: (p, S, Kc) => m.riemannDecode(p, S, Kc, 1, true), h: x => x < 600 ? 0.3 * P * 1100 : 0.075 * P * x, lo: 1, hi: 199000, n: 100, full: 4000 },
  135: { f: m.riemannClass, h: x => Math.max(0.125 * P * 1500, 0.045 * P * x), lo: 100, hi: 23000, n: 60, full: 4000, sweep: [100, 8100] },
  136: { f: m.riemannMod24, h: x => Math.max(0.033 * P * 4400, 0.008 * P * x), lo: 100, hi: 85000, n: 60, full: 4000, sweep: [4000, 85000] },
  137: { f: m.riemannMod840, h: x => 82, lo: 4000, hi: 79500, n: 60, sweep: [4000, 20000] },
  138: { f: m.riemannAuto, h: x => Math.max(0.033 * P * 4400, 0.045 * P * x), lo: 1, hi: 199000, n: 60, full: 4000 },
  139: { f: (p, S, Kc) => m.riemannDecode(p, S, Kc, 4, true), h: x => x < 600 ? 0.3 * P * 1100 : 0.045 * P * x, lo: 1, hi: 199000, n: 100, full: 4000 },
};
for (const k of [2, 4, 11, 25]) F['140k' + k] = { f: (p, S, Kc) => m.riemannDecode(p, S, Kc, k, true), h: x => x < 600 ? 0.3 * P * 1100 : 0.23 * rho(k) * P * x, lo: 600, hi: Math.min(300000, Math.floor(G[K - 1] / (0.23 * rho(k) * P)) - 900), n: 30, full: 3000 };
const run = (f, S, C, Kc) => { clear(); let bad = 0, tot = 0, first = ''; for (let p = Math.max(2, S); p <= S + C; p++) if (isP(p)) { tot++; let q = p + 1; while (!isP(q)) q++; const y = f(p, S, Kc); if (y !== q) { bad++; if (!first) first = `${p} -> ${y}`; } } return [bad, tot, first]; };
const want = (process.argv[2] || '134,135,136,137,138,139,140').split(',');
for (const name of Object.keys(F)) {
  if (!want.includes(name) && !(name.startsWith('140') && want.includes('140'))) continue;
  const c = F[name], t = Date.now(); let out = `${name}:`;
  if (c.full) { const r = run(c.f, 1, c.full, cnt(c.h(c.full + 500) * 1.02)); out += ` 1-${c.full}: ${r[0]} wrong of ${r[1]} primes ${r[2]};`; }
  if (c.sweep) { const [s0, s1] = c.sweep, r = run(c.f, s0, s1 - s0, cnt(c.h(s1 + 600) * 1.02)); out += ` ${s0}-${s1} in one run: ${r[0]} wrong of ${r[1]} ${r[2]};`; }
  let seed = 7 + name.length * 31 + Number(name.replace(/\D/g, '')); const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  let B = 0, T = 0; const w = [];
  for (let i = 0; i < c.n; i++) { const C = 150 + Math.floor(rnd() * 450), S = Math.floor(c.lo + rnd() * (c.hi - c.lo - C)); const r = run(c.f, S, C, cnt(c.h(S + C + 600) * 1.02)); B += r[0]; T += r[1]; if (r[0]) w.push(`S=${S} ${r[2]}`); }
  console.log(`${out} ${c.n} random circles ${c.lo}-${c.hi}: ${B} wrong of ${T} ${w.slice(0, 3).join('; ')} [${((Date.now() - t) / 1000).toFixed(0)} s]`);
}
