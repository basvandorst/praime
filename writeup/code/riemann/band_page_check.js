// Validates transformation 141 of prime.html (primes from a band of zeros, RIEMANN.md section 12).
// It runs the page's own code (extracted from prime.html), with Odlyzko's tables from writeup/data/riemann
// and the page's own zeros.
//   node band_page_check.js reach [bands] [ks]   decode from 3 with nothing known (8 histories, the page's setting) up to
//                                                15 band widths (the page's limit); the first wrong decision must equal the
//                                                number in the page's band table (a negative entry means clean to that number)
//   node band_page_check.js circles [bands]      the page's own mode: circles of 1,000 numbers with the primes below the
//                                                circle known; every next prime it returns is checked by trial division
// bands: comma list of slider positions 0-7 (default all); ks: comma list of k (default 1,2,3,4).
// PAGE=path overrides the page file (default: prime.html next to this script); OUT=dir writes each decided list of
// 'reach' to dir/band<i>_k<k>.txt (a line '# N=limit', then one prime per line), for the comparison with band_decoder.py in band_page_cross.sh.
const fs = require('fs'), path = require('path');
const html = fs.readFileSync(process.env.PAGE || path.join(__dirname, 'prime.html'), 'utf8');
const js = html.match(/<script>\n([\s\S]*)<\/script>/)[1];
const a = js.indexOf('const Z30='), b = js.indexOf('let sieveCache');
let BAND = '6';
const $ = id => ({ value: id === 'band' ? BAND : '1576', textContent: '' });
const m = new Function('$', js.slice(a, b) + ';return {BANDS,bandData,bandParse,bandSet,bandState,bandStep,bandFlush,riemannBand,bandStates,DEC_PRIMES,zetaZeros};')($);
const D = path.join(__dirname, '..', '..', 'data', 'riemann');
for (const f of ['odlyzko_zeros3.txt', 'odlyzko_zeros4.txt', 'odlyzko_zeros5.txt']) m.bandData.set(f, m.bandParse(fs.readFileSync(path.join(D, f), 'utf8')));
const isP = n => { if (n < 2) return false; if (n % 2 === 0) return n === 2; for (let q = 3; q * q <= n; q += 2) if (n % q === 0) return false; return true; };
const [mode, bl, kl] = process.argv.slice(2);
const bands = (bl || '0,1,2,3,4,5,6,7').split(',').map(Number), ks = (kl || '1,2,3,4').split(',').map(Number);
let bad = 0;
if (mode === 'reach') {
  for (const i of bands) {
    BAND = String(i);
    const t0 = Date.now(), bs = m.bandSet(i), W = bs.reduce((s, x) => s + x.width, 0), cap = Math.floor(15 * W);
    for (const k of ks) {
      const t = Date.now(), wheel = m.DEC_PRIMES.slice(0, k), st = m.bandState(bs, 3, wheel, 8, cap + 4000), sv = new Uint8Array(cap + 4001);
      for (let x = 2; x <= cap + 4000; x++) sv[x] = 1;
      for (let x = 2; x * x <= cap + 4000; x++) if (sv[x]) for (let y = x * x; y <= cap + 4000; y += x) sv[y] = 0;
      const errors = [], list = []; let count = 0, seen = 0;
      const scan = () => { let j = 0; for (const [n, v] of st.dec) { if (j++ < seen) continue; seen++; if (v) { count++; list.push(n); } if (!!v !== !!sv[n]) errors.push(n); if (errors.length >= 3) return true; } return false; };
      let stop = false;
      while (!stop && st.i < st.allc.length && st.allc[st.i] <= cap) { m.bandStep(st); stop = scan(); }
      if (!stop) { m.bandFlush(st); scan(); }
      if (process.env.OUT) fs.writeFileSync(path.join(process.env.OUT, `band${i}_k${k}.txt`), `# N=${cap}\n` + list.join('\n') + '\n');
      const got = errors.length ? errors[0] : -cap, want = m.BANDS[i].reach[k - 1], ok = got === want;
      if (!ok) bad++;
      console.log(`band ${i} (${m.BANDS[i].name}), width ${W.toFixed(0)}, k = ${k}: ${count} primes decided; first wrong [${errors.join(', ')}]; ` +
        `${errors.length ? 'reach ' + errors[0] + ' = ' + (errors[0] / W).toFixed(2) + ' x width' : 'clean to ' + cap} ; page says ${want}: ${ok ? 'ok' : 'MISMATCH'} [${((Date.now() - t) / 1000).toFixed(0)} s]`);
    }
  }
} else if (mode === 'circles') {
  let seed = 141; const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (const i of bands) {
    BAND = String(i);
    const bs = m.bandSet(i), W = bs.reduce((s, x) => s + x.width, 0), R = m.BANDS[i].reach.map(Math.abs);
    for (const k of ks) {
      const reach = R[k - 1], starts = [1, ...[0.25, 0.5, 0.75].map(f => Math.max(2, Math.floor(f * reach - 500 + 200 * rnd())))];
      const parts = [];
      for (const S of starts) {
        const t = Date.now(); m.bandStates.clear(); let wrong = 0, tot = 0, first = '';
        for (let p = Math.max(2, S); p < S + 1000 && p < reach - 400; p++) if (isP(p)) { tot++; let q = p + 1; while (!isP(q)) q++; const y = m.riemannBand(p, S, k); if (y !== q) { wrong++; if (!first) first = ` (first: ${p} -> ${y}, true ${q})`; } }
        if (wrong) bad++;
        parts.push(`S = ${S}: ${wrong} wrong of ${tot}${first} [${((Date.now() - t) / 1000).toFixed(1)} s]`);
      }
      console.log(`band ${i}, k = ${k}, circles of 1,000 below the reach ${reach}: ${parts.join('; ')}`);
    }
  }
}
console.log(bad ? `${bad} FAILED` : 'all checks passed');
