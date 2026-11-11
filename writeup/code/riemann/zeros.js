// Computes the first K zeta zeros with the zero finder of prime.html
// (Riemann-Siegel with C0, C1, C2; sign scan with close-pair check; Taylor step)
// and writes them as float64 to a file.
// usage: node zeros.js K out.f64
const fs = require('fs'), path = require('path');
const html = fs.readFileSync(path.join(__dirname, 'prime.html'), 'utf8');
const js = html.match(/<script>\n([\s\S]*)<\/script>/)[1];
const a = js.indexOf('const Z30='), b = js.indexOf('let sieveCache');
const page = new Function('$', js.slice(a, b) + ';return {zetaZeros};')(() => ({ value: 1576 }));
const K = Number(process.argv[2] || 10142), out = process.argv[3];
if (!out) { console.log('usage: node zeros.js K out.f64'); process.exit(1); }
const t = Date.now(), G = page.zetaZeros(K).slice(0, K);
fs.writeFileSync(out, Buffer.from(Float64Array.from(G).buffer));
console.log(`${G.length} zeros, the last at ${G[G.length - 1].toFixed(6)} [${((Date.now() - t) / 1000).toFixed(1)} s]`);
