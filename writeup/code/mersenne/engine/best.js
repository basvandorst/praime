'use strict';
/* best.js: p -> the next Mersenne exponent after p (the next prime q > p with 2^q - 1 prime).
   Plain functions, no libraries, no Node APIs in the core: load as <script> (global MersenneNext), importScripts()
   in a Web Worker, or require() in Node.  nextMersenneExponent(p[, opts]), testExponent(q[, opts]).

   Per prime q:
   1. trial factoring: f = 2kq + 1, k <= K(q) = 1e5 (q/7000)^2, only k with f = +-1 mod 8, k sieved by the odd
      primes below 2^14; test 2^q = 1 mod f in double arithmetic (f < 2^38, products split into 13- or 16-bit
      chunks so that every product and every reduction is exact). A factor found proves 2^q - 1 composite.
   2. Lucas-Lehmer (theorem): s = 4, s <- s^2 - 2 mod 2^q - 1, q - 2 times; 2^q - 1 is prime iff s = 0.
      q < 4200: BigInt square and fold (exact). q >= 4200: IBDWT (Crandall-Fagin) squaring in Float64Array:
      n = 2^m or 3*2^m digits of floor/ceil(q/n) bits, weights 2^(ceil(qi/n) - qi/n), balanced digits, real-signal
      FFT through a complex radix-4 (+ one radix-3 pass) FFT of length n/2, at most 25 - log2(n)/2 bits per digit.
      The largest rounding error of every iteration is checked; above 0.35 the test restarts at twice the length.
      A zero residue from the FFT is re-checked with the exact BigInt test (opts.confirm, default on), so every
      exponent returned rests on exact arithmetic; a composite verdict from the FFT rests on the measured rounding
      error (max ~0.03 over whole tests, fatal at 0.5), not on a proof. */

const RND = 6755399441055744; // 1.5 * 2^52: (x + RND) - RND rounds x to the nearest integer for |x| < 2^51

function isPrimeSmall(n) {
  if (n < 2) return false;
  if (n % 2 === 0) return n === 2;
  if (n % 3 === 0) return n === 3;
  for (let d = 5; d * d <= n; d += 6) if (n % d === 0 || n % (d + 2) === 0) return false;
  return true;
}

/* ---------------- Lucas-Lehmer with BigInt ---------------- */
function llResidueBigInt(q) { // returns s_{q-2} mod 2^q - 1 as a BigInt in [0, M)
  if (q === 2) return 0n;
  const Q = BigInt(q), M = (1n << Q) - 1n;
  let s = 4n;
  for (let i = 0; i < q - 2; i++) {
    s = s * s - 2n;
    s = (s & M) + (s >> Q);
    if (s >= M) s -= M;
  }
  return s;
}

/* ---------------- IBDWT squaring mod 2^q - 1 ---------------- */
// FFT lengths n = 2^k and 3*2^k. Largest average bits per digit for length n, measured: the max rounding error over
// 2000 LL iterations is 0.02-0.05 at this size (n = 128..4096) and doubles with every further half bit (0.5 is fatal).
function maxBitsPerDigit(n) { return 25 - 0.5 * Math.log2(n); }
function fftLengthFor(q) {
  for (let n = 64; ; n *= 2) {
    if (q / n <= maxBitsPerDigit(n)) return n;
    if (n >= 256 && q / (1.5 * n) <= maxBitsPerDigit(1.5 * n)) return 1.5 * n;
  }
}

function makePlan(q, n) {
  const N = n >> 1, R3 = N % 3 === 0, M = R3 ? N / 3 : N, LOG = Math.round(Math.log2(M));
  const c = new Float64Array(n + 1);                       // c_i = ceil(q i / n): bit position of digit i
  for (let i = 0; i <= n; i++) c[i] = Math.floor((q * i + n - 1) / n);
  const wt = new Float64Array(n), iwt = new Float64Array(n), base = new Float64Array(n), ibase = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const fr = (c[i] * n - q * i) / n;                       // in [0, 1)
    wt[i] = Math.pow(2, fr);
    iwt[i] = 1 / (wt[i] * N);                                // unweight and the 1/N of the inverse FFT
    const b = c[i + 1] - c[i];
    base[i] = Math.pow(2, b); ibase[i] = Math.pow(2, -b);
  }
  const cs = new Float64Array(N), sn = new Float64Array(N);   // cos, sin(2 pi t/N), full circle
  for (let t = 0; t < N; t++) { cs[t] = Math.cos(2 * Math.PI * t / N); sn[t] = Math.sin(2 * Math.PI * t / N); }
  // pos[k] = where Z_k sits after the forward transform: (k mod 3) M + bitrev_M(floor(k/3)) (radix 3), else bitrev(k)
  const br = new Int32Array(N);
  for (let k = 0; k < N; k++) {
    const k1 = R3 ? k % 3 : 0, k2 = R3 ? (k - k1) / 3 : k;
    let r = 0; for (let b = 0; b < LOG; b++) if (k2 >> b & 1) r |= 1 << (LOG - 1 - b);
    br[k] = k1 * M + r;
  }
  return { q, n, N, M, R3, c, wt, iwt, base, ibase, cs, sn, br, re: new Float64Array(N), im: new Float64Array(N) };
}

// Complex FFTs of length N (power of 2, >= 4) on re/im[off .. off+N), twiddle e^{2 pi i t/N} = table index t*ts.
// Radix-4 passes = two fused radix-2 passes, so the output order (bit-reversed) is the same as radix 2.
// DIF pair (h, h/2), block 2h, j < h/2, w = e^{-2 pi i j/(2h)}:
//   y0 = (x0+x2)+(x1+x3), y1 = w^2 [(x0+x2)-(x1+x3)], y2 = w [(x0-x2) - i(x1-x3)], y3 = w^3 [(x0-x2) + i(x1-x3)]
function fftDIF4(re, im, off, N, cs, sn, ts) {
  let h = N >> 1, st = ts;
  const end = off + N;
  if ((31 - Math.clz32(N)) & 1) {                             // odd log2 N: one radix-2 pass first
    for (let j = 0; j < h; j++) {
      const a = off + j, b = a + h, wr = cs[j * ts], wi = -sn[j * ts];
      const xr = re[a] - re[b], xi = im[a] - im[b];
      re[a] += re[b]; im[a] += im[b];
      re[b] = xr * wr - xi * wi; im[b] = xr * wi + xi * wr;
    }
    h >>= 1; st <<= 1;
  }
  for (; h >= 4; h >>= 2, st <<= 2) {
    const q4 = h >> 1;
    for (let s = off; s < end; s += 2 * h) {
      for (let j = 0, t = 0; j < q4; j++, t += st) {
        const a0 = s + j, a1 = a0 + q4, a2 = a1 + q4, a3 = a2 + q4;
        const x0r = re[a0], x0i = im[a0], x1r = re[a1], x1i = im[a1], x2r = re[a2], x2i = im[a2], x3r = re[a3], x3i = im[a3];
        const sr = x0r + x2r, si = x0i + x2i, dr = x0r - x2r, di = x0i - x2i;
        const tr = x1r + x3r, ti = x1i + x3i, er = x1r - x3r, ei = x1i - x3i;
        re[a0] = sr + tr; im[a0] = si + ti;
        const ur = sr - tr, ui = si - ti;                     // * w^2
        const w2r = cs[2 * t], w2i = -sn[2 * t];
        re[a1] = ur * w2r - ui * w2i; im[a1] = ur * w2i + ui * w2r;
        const vr = dr + ei, vi = di - er;                     // (d - i e) * w
        const w1r = cs[t], w1i = -sn[t];
        re[a2] = vr * w1r - vi * w1i; im[a2] = vr * w1i + vi * w1r;
        const zr = dr - ei, zi = di + er;                     // (d + i e) * w^3
        const w3r = cs[3 * t], w3i = -sn[3 * t];
        re[a3] = zr * w3r - zi * w3i; im[a3] = zr * w3i + zi * w3r;
      }
    }
  }
  if (h === 2) {                                              // last pair (2, 1): all twiddles 1
    for (let a0 = off; a0 < end; a0 += 4) {
      const a1 = a0 + 1, a2 = a0 + 2, a3 = a0 + 3;
      const x0r = re[a0], x0i = im[a0], x1r = re[a1], x1i = im[a1], x2r = re[a2], x2i = im[a2], x3r = re[a3], x3i = im[a3];
      const sr = x0r + x2r, si = x0i + x2i, dr = x0r - x2r, di = x0i - x2i;
      const tr = x1r + x3r, ti = x1i + x3i, er = x1r - x3r, ei = x1i - x3i;
      re[a0] = sr + tr; im[a0] = si + ti; re[a1] = sr - tr; im[a1] = si - ti;
      re[a2] = dr + ei; im[a2] = di - er; re[a3] = dr - ei; im[a3] = di + er;
    }
  }
}
// DIT pair (h, 2h), block 4h, j < h, u = e^{+2 pi i j/(4h)}: B0,1 = x0 +- u^2 x1, C0,1 = u x2 +- u^3 x3,
//   y0 = B0 + C0, y2 = B0 - C0, y1 = B1 + i C1, y3 = B1 - i C1
function fftDIT4(re, im, off, N, cs, sn, ts) {
  const LOG = 31 - Math.clz32(N), end = off + N;
  for (let a0 = off; a0 < end; a0 += 4) {                     // first pair (1, 2): all twiddles 1
    const a1 = a0 + 1, a2 = a0 + 2, a3 = a0 + 3;
    const x0r = re[a0], x0i = im[a0], x1r = re[a1], x1i = im[a1], x2r = re[a2], x2i = im[a2], x3r = re[a3], x3i = im[a3];
    const b0r = x0r + x1r, b0i = x0i + x1i, b1r = x0r - x1r, b1i = x0i - x1i;
    const c0r = x2r + x3r, c0i = x2i + x3i, c1r = x2r - x3r, c1i = x2i - x3i;
    re[a0] = b0r + c0r; im[a0] = b0i + c0i; re[a2] = b0r - c0r; im[a2] = b0i - c0i;
    re[a1] = b1r - c1i; im[a1] = b1i + c1r; re[a3] = b1r + c1i; im[a3] = b1i - c1r;
  }
  for (let h = 4, st = (N >> 4) * ts; 4 * h <= N; h <<= 2, st >>= 2) {
    for (let s = off; s < end; s += 4 * h) {
      for (let j = 0, t = 0; j < h; j++, t += st) {
        const a0 = s + j, a1 = a0 + h, a2 = a1 + h, a3 = a2 + h;
        const x0r = re[a0], x0i = im[a0];
        const w2r = cs[2 * t], w2i = sn[2 * t], w1r = cs[t], w1i = sn[t], w3r = cs[3 * t], w3i = sn[3 * t];
        const y1r = re[a1] * w2r - im[a1] * w2i, y1i = re[a1] * w2i + im[a1] * w2r;
        const y2r = re[a2] * w1r - im[a2] * w1i, y2i = re[a2] * w1i + im[a2] * w1r;
        const y3r = re[a3] * w3r - im[a3] * w3i, y3i = re[a3] * w3i + im[a3] * w3r;
        const b0r = x0r + y1r, b0i = x0i + y1i, b1r = x0r - y1r, b1i = x0i - y1i;
        const c0r = y2r + y3r, c0i = y2i + y3i, c1r = y2r - y3r, c1i = y2i - y3i;
        re[a0] = b0r + c0r; im[a0] = b0i + c0i; re[a2] = b0r - c0r; im[a2] = b0i - c0i;
        re[a1] = b1r - c1i; im[a1] = b1i + c1r; re[a3] = b1r + c1i; im[a3] = b1i - c1r;
      }
    }
  }
  if (LOG & 1) {                                              // odd log2 N: last radix-2 pass, half-span N/2
    const H = N >> 1;
    for (let j = 0; j < H; j++) {
      const a = off + j, b = a + H, wr = cs[j * ts], wi = sn[j * ts];
      const tr = re[b] * wr - im[b] * wi, ti = re[b] * wi + im[b] * wr;
      re[b] = re[a] - tr; im[b] = im[a] - ti;
      re[a] += tr; im[a] += ti;
    }
  }
}
const S3 = Math.sqrt(3) / 2;
// Forward transform of length N = M or 3M. Radix 3 (DIF): Y0 = x0+x1+x2, Y1,2 = x0 - (x1+x2)/2 -+ i S3 (x1-x2),
// then Y_k1[j] *= e^{-2 pi i j k1/N}, and each third gets a length-M FFT.
function fwd(P) {
  const { re, im, N, M, cs, sn } = P;
  if (!P.R3) { fftDIF4(re, im, 0, N, cs, sn, 1); return; }
  for (let j = 0; j < M; j++) {
    const b = j + M, c = b + M;
    const tr = re[b] + re[c], ti = im[b] + im[c], dr = re[b] - re[c], di = im[b] - im[c];
    const mr = re[j] - 0.5 * tr, mi = im[j] - 0.5 * ti;
    re[j] += tr; im[j] += ti;
    const y1r = mr + S3 * di, y1i = mi - S3 * dr, y2r = mr - S3 * di, y2i = mi + S3 * dr;
    const w1r = cs[j], w1i = -sn[j], w2r = cs[2 * j], w2i = -sn[2 * j];
    re[b] = y1r * w1r - y1i * w1i; im[b] = y1r * w1i + y1i * w1r;
    re[c] = y2r * w2r - y2i * w2i; im[c] = y2r * w2i + y2i * w2r;
  }
  fftDIF4(re, im, 0, M, cs, sn, 3); fftDIF4(re, im, M, M, cs, sn, 3); fftDIF4(re, im, 2 * M, M, cs, sn, 3);
}
// Inverse (unscaled): length-M inverse FFTs on the thirds, then v_k1 = u_k1[j] e^{+2 pi i j k1/N},
// z0 = u0+v1+v2, z1,2 = u0 - (v1+v2)/2 +- i S3 (v1-v2)
function inv(P) {
  const { re, im, N, M, cs, sn } = P;
  if (!P.R3) { fftDIT4(re, im, 0, N, cs, sn, 1); return; }
  fftDIT4(re, im, 0, M, cs, sn, 3); fftDIT4(re, im, M, M, cs, sn, 3); fftDIT4(re, im, 2 * M, M, cs, sn, 3);
  for (let j = 0; j < M; j++) {
    const b = j + M, c = b + M;
    const w1r = cs[j], w1i = sn[j], w2r = cs[2 * j], w2i = sn[2 * j];
    const v1r = re[b] * w1r - im[b] * w1i, v1i = re[b] * w1i + im[b] * w1r;
    const v2r = re[c] * w2r - im[c] * w2i, v2i = re[c] * w2i + im[c] * w2r;
    const tr = v1r + v2r, ti = v1i + v2i, dr = v1r - v2r, di = v1i - v2i;
    const mr = re[j] - 0.5 * tr, mi = im[j] - 0.5 * ti;
    re[j] += tr; im[j] += ti;
    re[b] = mr - S3 * di; im[b] = mi + S3 * dr;
    re[c] = mr + S3 * di; im[c] = mi - S3 * dr;
  }
}

/* x <- mul * x^2 - sub (mod 2^q - 1), digits x[0..n-1] balanced. Returns the largest rounding error. */
function sqrStep(P, x, sub, mul) {
  const { n, N, wt, iwt, base, ibase, cs, sn, br, re, im } = P;
  for (let j = 0; j < N; j++) { re[j] = x[2 * j] * wt[2 * j]; im[j] = x[2 * j + 1] * wt[2 * j + 1]; }
  fwd(P);
  // real-signal squaring: Z_k = E_k + i O_k (E, O = DFTs of the even and odd digits), Y_k^2 folded back into
  // Z'_k = E_k^2 + w^k O_k^2 + 2i E_k O_k with w = e^{-2 pi i/N}; k and N-k are processed together.
  const H = N >> 1;
  for (let k = 0; k <= H; k++) {
    const m = k === 0 ? 0 : N - k, pk = br[k], pm = br[m];
    const ar = re[pk], ai = im[pk], b_r = re[pm], b_i = im[pm];
    const er = 0.5 * (ar + b_r), ei = 0.5 * (ai - b_i);       // E_k
    const or = 0.5 * (ai + b_i), oi = 0.5 * (b_r - ar);       // O_k
    const wr = cs[k], wi = -sn[k];
    const e2r = er * er - ei * ei, e2i = 2 * er * ei;
    const o2r = or * or - oi * oi, o2i = 2 * or * oi;
    const pr = e2r + wr * o2r - wi * o2i, pi = e2i + wr * o2i + wi * o2r;  // P = E^2 + w O^2
    const qr = er * or - ei * oi, qi = er * oi + ei * or;                 // Q = E O
    re[pk] = pr - 2 * qi; im[pk] = pi + 2 * qr;
    if (m !== k) { re[pm] = pr + 2 * qi; im[pm] = -pi + 2 * qr; }
  }
  inv(P);
  let carry = -sub, err = 0;
  for (let j = 0; j < N; j++) {
    let i = 2 * j;
    let v = re[j] * iwt[i], r = (v + RND) - RND, e = v - r;
    if (e < 0) e = -e; if (e > err) err = e;
    let t = r * mul + carry, hi = (t * ibase[i] + RND) - RND;
    x[i] = t - hi * base[i]; carry = hi;
    i++;
    v = im[j] * iwt[i]; r = (v + RND) - RND; e = v - r;
    if (e < 0) e = -e; if (e > err) err = e;
    t = r * mul + carry; hi = (t * ibase[i] + RND) - RND;
    x[i] = t - hi * base[i]; carry = hi;
  }
  for (let i = 0; carry !== 0 && i < n; i++) {               // 2^q = 1: the top carry wraps to digit 0
    const t = x[i] + carry, hi = (t * ibase[i] + RND) - RND;
    x[i] = t - hi * base[i]; carry = hi;
  }
  return carry === 0 ? err : Infinity;                        // NaN or a runaway carry is reported as a failure
}

function digitsToBigInt(P, x) { // value mod 2^q - 1 in [0, M)
  const M = (1n << BigInt(P.q)) - 1n;
  let v = 0n;
  for (let i = P.n - 1; i >= 0; i--) v += BigInt(x[i]) << BigInt(P.c[i]);
  v %= M; if (v < 0n) v += M;
  return v;
}
function bigIntToDigits(P, v) {
  const x = new Float64Array(P.n);
  for (let i = 0; i < P.n; i++) x[i] = Number((v >> BigInt(P.c[i])) & ((1n << BigInt(P.c[i + 1] - P.c[i])) - 1n));
  return x;
}

const ERR_MAX = 0.35;
const planCache = new Map();
function getPlan(q, n) { const key = q + ':' + n; let P = planCache.get(key); if (!P) { if (planCache.size > 8) planCache.clear(); P = makePlan(q, n); planCache.set(key, P); } return P; }

// Lucas-Lehmer by IBDWT. Returns { residue (BigInt), n, maxErr }.
function llResidueFFT(q, n) {
  n = n || fftLengthFor(q);
  for (;;) {
    const P = getPlan(q, n), x = new Float64Array(n);
    x[0] = 4;
    let maxErr = 0, ok = true;
    for (let it = 0; it < q - 2; it++) {
      const e = sqrStep(P, x, 2, 1);
      if (!(e <= ERR_MAX)) { ok = false; break; }             // also catches NaN
      if (e > maxErr) maxErr = e;
    }
    if (ok) return { residue: digitsToBigInt(P, x), n, maxErr };
    n *= 2;                                                   // rounding too close to 1/2: redo at twice the length
  }
}

/* ---------------- P-1 stage 1 ---------------- */
function primesUpTo(B) { const s = new Uint8Array(B + 1), out = []; for (let i = 2; i <= B; i++) if (!s[i]) { out.push(i); for (let j = i * i; j <= B; j += i) s[j] = 1; } return out; }
function gcdBig(a, b) { while (b) { const t = a % b; a = b; b = t; } return a; }
// returns a proper factor of 2^q - 1 found by P-1 stage 1 with bound B1, else 0n. Exponent 2 q E, E = prod p^floor(log_p B1).
function pm1Stage1(q, B1, useFFT) {
  const M = (1n << BigInt(q)) - 1n;
  let E = 2n * BigInt(q);
  for (const p of primesUpTo(B1)) { let pk = p; while (pk * p <= B1) pk *= p; E *= BigInt(pk); }
  const bits = E.toString(2);
  let x;
  if (useFFT) {
    const n = fftLengthFor(q), P = getPlan(q, n), d = new Float64Array(n); d[0] = 3;
    for (let i = 1; i < bits.length; i++) sqrStep(P, d, 0, bits.charCodeAt(i) === 49 ? 3 : 1);
    x = digitsToBigInt(P, d);
  } else {
    x = 3n; const Q = BigInt(q);
    for (let i = 1; i < bits.length; i++) {
      x = x * x; if (bits.charCodeAt(i) === 49) x *= 3n;
      x = (x & M) + (x >> Q); if (x >= M) x -= M;
    }
  }
  const g = gcdBig(M, (x - 1n + M) % M);
  return (g > 1n && g < M) ? g : 0n;
}

/* ---------------- trial factoring, f = 2kq + 1 < 2^38 ---------------- */
// 2^q mod f. The top bits of q are preloaded (2^v, v < 32). Squaring r*r mod f splits r into chunks so every product
// is exact: f < 2^32: two 16-bit chunks (x < 2^49); f < 2^38: three 13-bit chunks (x = acc*2^13 + r*chunk < 2^52).
// Reduction x - floor(x/f)*f is exact (x < 2^52), and floor is off by at most 1 (quotient < 2^14), fixed by one step.
function pow2mod(q, f) {
  const finv = 1 / f;
  let bit = 31 - Math.clz32(q), v = 0;
  while (bit >= 0 && 2 * v + ((q >>> bit) & 1) < 32) { v = 2 * v + ((q >>> bit) & 1); bit--; }
  let r = Math.pow(2, v) % f;
  if (f < 4294967296) {
    for (; bit >= 0; bit--) {
      const rh = Math.floor(r * 1.52587890625e-5), rl = r - rh * 65536;
      let x = r * rh, t = Math.floor(x * finv); x -= t * f; if (x < 0) x += f; else if (x >= f) x -= f;
      x = x * 65536 + r * rl; t = Math.floor(x * finv); x -= t * f; if (x < 0) x += f; else if (x >= f) x -= f;
      r = x;
      if ((q >>> bit) & 1) { r += r; if (r >= f) r -= f; }
    }
  } else {
    for (; bit >= 0; bit--) {
      const r2 = Math.floor(r * 1.4901161193847656e-8), rr = r - r2 * 67108864;   // r = r2 2^26 + r1 2^13 + r0
      const r1 = Math.floor(rr * 1.220703125e-4), r0 = rr - r1 * 8192;
      let x = r * r2, t = Math.floor(x * finv); x -= t * f; if (x < 0) x += f; else if (x >= f) x -= f;
      x = x * 8192 + r * r1; t = Math.floor(x * finv); x -= t * f; if (x < 0) x += f; else if (x >= f) x -= f;
      x = x * 8192 + r * r0; t = Math.floor(x * finv); x -= t * f; if (x < 0) x += f; else if (x >= f) x -= f;
      r = x;
      if ((q >>> bit) & 1) { r += r; if (r >= f) r -= f; }
    }
  }
  return r;
}
function pow2modBig(q, f) { let r = 1n, b = 2n, e = q; const F = BigInt(f); while (e > 0) { if (e & 1) r = r * b % F; b = b * b % F; e >>>= 1; } return Number(r); }

let SIEVE_PRIMES = null;
function sievePrimes() { if (!SIEVE_PRIMES) SIEVE_PRIMES = primesUpTo(1 << 14).slice(1); return SIEVE_PRIMES; } // odd primes < 2^14
function invMod(a, m) { let t = 0, nt = 1, r = m, nr = a % m; while (nr) { const qq = Math.floor(r / nr); let x = t - qq * nt; t = nt; nt = x; x = r - qq * nr; r = nr; nr = x; } return t < 0 ? t + m : t; }

// A factor f = 2kq + 1 (k <= K) of 2^q - 1, else 0. Only k in the two classes mod 4 with f = +-1 mod 8 (k = 4t + c),
// sieved by the odd primes below 2^14 (prime s kills t = ((-(2q)^-1 - c) 4^-1) mod s; f = s itself is kept).
// Segments of 8192 t, both classes per segment, so small factors are met first.
function trialFactor(q, K) {
  if (K < 1) return 0;
  const Mq = q < 53 ? Math.pow(2, q) - 1 : Infinity;
  const classes = (q & 3) === 1 ? [0, 3] : [0, 1];
  const SP = sievePrimes(), NS = SP.length, SEG = 8192, seg = new Uint8Array(SEG);
  const tmax = Math.floor((K - Math.min(...classes)) / 4) + 1;
  const roots = classes.map(() => new Int32Array(NS)), self = classes.map(() => new Float64Array(NS).fill(-1));
  classes.forEach((c, ci) => {
    for (let i = 0; i < NS; i++) {
      const s = SP[i];
      if (s === q) { roots[ci][i] = -1; continue; }
      const ks = (s - invMod((2 * q) % s, s)) % s;
      roots[ci][i] = (((ks - c) % s + s) % s) * invMod(4, s) % s;
      if (s > 2 * q && (s - 1) % (2 * q) === 0) { const k = (s - 1) / (2 * q); if (k % 4 === c) self[ci][i] = (k - c) / 4; }
    }
  });
  for (let T0 = 0; T0 < tmax; T0 += SEG) {
    const len = Math.min(SEG, tmax - T0);
    for (let ci = 0; ci < classes.length; ci++) {
      const c = classes[ci], R = roots[ci], S = self[ci];
      seg.fill(0);
      for (let i = 0; i < NS; i++) {
        const s = SP[i], r0 = R[i]; if (r0 < 0) continue;
        let t = r0 - (T0 % s); if (t < 0) t += s;
        for (; t < len; t += s) seg[t] = 1;
        if (S[i] >= T0 && S[i] < T0 + len) seg[S[i] - T0] = 0;     // f = s is a prime candidate itself
      }
      for (let t = 0; t < len; t++) {
        if (seg[t]) continue;
        const k = 4 * (T0 + t) + c;
        if (k < 1) continue;
        if (k > K) break;
        const f = 2 * k * q + 1;
        if (f >= Mq) break;
        const r = f < 274877906944 ? pow2mod(q, f) : pow2modBig(q, f);
        if (r === 1) return f;
      }
    }
  }
  return 0;
}

/* ---------------- the transformation ---------------- */
// Defaults (measured in Node 22, see the report): FFT from q = 4200 (BigInt is faster below); trial-factoring depth
// K(q) = 1e5 (q/7000)^2, capped so that f < 2^38 (the TF optimum moves one bit deeper when the LL cost doubles, and LL
// cost grows like q^2); P-1 off (B1 = 0: a net loss at these sizes); confirm = true re-runs LL in exact BigInt
// arithmetic whenever the FFT returns residue 0, so every reported exponent rests on exact arithmetic only.
const DEFAULTS = { K: null, FFT_MIN: 4200, B1: 0, confirm: true };
function tfDepth(q) { return Math.max(1000, Math.min(Math.round(1e5 * (q / 7000) ** 2), Math.floor(2.7e11 / (2 * q)))); }
// Is 2^q - 1 prime (q prime)? Returns { prime, how, factor?, n?, maxErr? }
function testExponent(q, opts) {
  const o = Object.assign({}, DEFAULTS, opts);
  if (q === 2) return { prime: true, how: 'q=2' };
  const f = trialFactor(q, o.K == null ? tfDepth(q) : o.K);
  if (f) return { prime: false, how: 'tf', factor: f };
  if (o.B1 > 0 && q >= o.FFT_MIN) { const g = pm1Stage1(q, o.B1, true); if (g) return { prime: false, how: 'pm1', factor: g }; }
  if (q < o.FFT_MIN) return { prime: llResidueBigInt(q) === 0n, how: 'll-bigint' };
  const r = llResidueFFT(q);
  let prime = r.residue === 0n, how = 'll-fft';
  if (prime && o.confirm) { prime = llResidueBigInt(q) === 0n; how = 'll-fft+bigint'; }
  return { prime, how, n: r.n, maxErr: r.maxErr };
}
function nextMersenneExponent(p, opts) {
  for (let q = Math.max(2, Math.floor(p) + 1); ; q++) {
    if (!isPrimeSmall(q)) continue;
    if (testExponent(q, opts).prime) return q;
  }
}

const api = { nextMersenneExponent, testExponent, isPrimeSmall, llResidueBigInt, llResidueFFT, makePlan, sqrStep,
  digitsToBigInt, bigIntToDigits, fftLengthFor, maxBitsPerDigit, trialFactor, tfDepth, pow2mod, pow2modBig, pm1Stage1, gcdBig, DEFAULTS };
if (typeof module !== 'undefined' && module.exports) module.exports = api; else globalThis.MersenneNext = api;
