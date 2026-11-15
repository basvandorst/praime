// nextExponentForecast(p): heuristic forecast of the first Mersenne exponent above p.
//
// Model (Wagstaff heuristic, a Poisson process over the primes): a prime q is an exponent with
// rate w(q) = e^gamma * log2(a q) / q, a = 6 if q = 1 (mod 4), a = 2 if q = 3 (mod 4), and
//     P(next > x) = exp(-sum_{p < q <= x, q prime} w(q)).
// Quantile at level l: the smallest x with sum_{p < q <= x} w(q) >= -ln(1 - l).
//
// p < 1e6 : the exact sum over the primes (small sieve, cached); the bounds are primes.
// p >= 1e6: the closed form of the integral (prime density 1/ln q, half the primes in each class)
//     Lambda(p, y) = e^gamma * [ log2(y/p) + log2(sqrt 12) * ln(ln y / ln p) ],
//   inverted by Newton's method. Checked against the exact prime sums for p in [1e6, 2.9e8]:
//   relative error below 1e-3 (see js_compare.txt).
//
// Memoryless: the forecast given that no exponent lies in (p, T] is nextExponentForecast(T).
// Returns { median, lo50, hi50, lo90, hi90, exact } (50% = quartiles, 90% = 5% and 95% points).
(function (root) {
  'use strict';
  const EG = Math.exp(0.5772156649015329);   // e^gamma
  const C12 = Math.log2(Math.sqrt(12));      // log2(sqrt 12): mean of log2 6 and log2 2
  const LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95];
  const TARGETS = LEVELS.map(l => -Math.log1p(-l));
  const EXACT_BELOW = 1e6;
  let comp = null, compN = 0;                // cached sieve: comp[i] = 1 if i is not prime

  function sieveTo(n) {
    if (n <= compN) return comp;
    const s = new Uint8Array(n + 1);
    s[0] = 1; s[1] = 1;
    for (let i = 2; i * i <= n; i++) if (!s[i]) for (let j = i * i; j <= n; j += i) s[j] = 1;
    comp = s; compN = n;
    return s;
  }

  function exactBounds(p) {
    for (let n = Math.ceil(4 * p) + 200; ; n *= 2) {
      const s = sieveTo(n), out = [];
      let S = 0, k = 0;
      for (let q = Math.floor(p) + 1; q <= n && k < TARGETS.length; q++) {
        if (s[q]) continue;
        S += EG * Math.log2(((q & 3) === 1 ? 6 : 2) * q) / q;
        while (k < TARGETS.length && S >= TARGETS[k]) { out.push(q); k++; }
      }
      if (k === TARGETS.length) return out;
    }
  }

  function closedBounds(p) {
    const lp = Math.log(p);
    return TARGETS.map(L => {
      let y = p * Math.pow(2, L / EG);       // start: the pure log-Poisson answer
      for (let i = 0; i < 50; i++) {
        const f = EG * (Math.log2(y / p) + C12 * Math.log(Math.log(y) / lp)) - L;
        const dfdy = EG / Math.LN2 * (1 + C12 * Math.LN2 / Math.log(y)) / y;
        const dy = f / dfdy;
        y -= dy;
        if (Math.abs(dy) < 1e-12 * y) break;
      }
      return Math.round(y);
    });
  }

  function nextExponentForecast(p) {
    if (!(p >= 2)) throw new RangeError('nextExponentForecast: p must be >= 2');
    const exact = p < EXACT_BELOW;
    const b = exact ? exactBounds(p) : closedBounds(p);
    return { median: b[2], lo50: b[1], hi50: b[3], lo90: b[0], hi90: b[4], exact };
  }

  root.nextExponentForecast = nextExponentForecast;
  if (typeof module !== 'undefined' && module.exports) module.exports = { nextExponentForecast };
})(typeof globalThis !== 'undefined' ? globalThis : this);
