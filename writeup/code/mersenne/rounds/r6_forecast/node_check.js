
const { nextExponentForecast } = require(process.argv[2]);
const pts = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const t0 = Date.now();
const out = pts.map(p => { const r = nextExponentForecast(p); return [r.lo90, r.lo50, r.median, r.hi50, r.hi90, r.exact]; });
console.log(JSON.stringify({ out, ms: Date.now() - t0 }));
