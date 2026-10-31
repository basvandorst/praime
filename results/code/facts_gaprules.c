// Backtest of 100 next-prime gap rules (20 base rules x 5 sight levels).
// Train on primes in [X-T, X), test on the first N primes >= X.
// usage: gapbt X T N   -> CSV: rule,y,hits,short,long,runmean,runmax
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;

static u64 mulmod(u64 a, u64 b, u64 m) { return (u128)a * b % m; }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static int isprime(u64 n) {
  if (n < 2) return 0;
  static const u64 sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (int i = 0; i < 12; i++) { if (n == sp[i]) return 1; if (n % sp[i] == 0) return 0; }
  u64 d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  for (int i = 0; i < 12; i++) {
    u64 x = powmod(sp[i], d, n); if (x == 1 || x == n - 1) continue;
    int ok = 0; for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
    if (!ok) return 0;
  }
  return 1;
}
static int SMALL[200], nsmall;   // primes <= 1000
static int hasf(u64 c, int y) { for (int i = 0; i < nsmall && SMALL[i] <= y; i++) { u64 q = SMALL[i]; if (q * q > c) return 0; if (c % q == 0) return 1; } return 0; }
static const int WP[] = {2, 3, 5, 7, 11, 13};
static int wheelgap(u64 n, int k) { int g = 1; for (;; g++) { u64 m = n + g; int bad = 0; for (int i = 0; i < k; i++) if (m % WP[i] == 0 && m != (u64)WP[i]) { bad = 1; break; } if (!bad) return g; } }

#define NR 20
static const int YS[5] = {0, 7, 13, 100, 1000};
static int *tab30, *tab210, *tab2310, *tab210g;   // learned most-frequent gap
static int evenr(double v) { int g = 2 * (int)lround(v / 2); return g < 2 ? 2 : g; }
static int cmpi(const void *a, const void *b) { return *(int *)a - *(int *)b; }

// base gap of rule r at prime p with gap history h[0..nh-1] (h[nh-1] = last gap)
static int base(int r, u64 p, const int *h, int nh, int *wheelk) {
  *wheelk = 0;
  int g1 = nh ? h[nh - 1] : 2;
  switch (r) {
    case 0: return 2;
    case 1: return 6;
    case 2: *wheelk = 2; return wheelgap(p, 2);
    case 3: *wheelk = 3; return wheelgap(p, 3);
    case 4: *wheelk = 4; return wheelgap(p, 4);
    case 5: *wheelk = 5; return wheelgap(p, 5);
    case 6: *wheelk = 6; return wheelgap(p, 6);
    case 7: return g1;
    case 8: return g1 + 2;
    case 9: return nh >= 2 ? h[nh - 2] : 2;
    case 10: { int m = g1; for (int i = 1; i <= 3 && i <= nh; i++) if (h[nh - i] < m) m = h[nh - i]; return m; }
    case 11: { int k = nh < 20 ? nh : 20; double s = 0; for (int i = 1; i <= k; i++) s += h[nh - i]; return k ? evenr(s / k) : 2; }
    case 12: { int k = nh < 20 ? nh : 20, t[20]; if (!k) return 2; for (int i = 0; i < k; i++) t[i] = h[nh - 1 - i]; qsort(t, k, sizeof(int), cmpi); return evenr(t[k / 2]); }
    case 13: { int k = nh < 20 ? nh : 20, best = 2, bc = 0; for (int i = 1; i <= k; i++) { int c = 0; for (int j = 1; j <= k; j++) if (h[nh - j] == h[nh - i]) c++; if (c > bc) { bc = c; best = h[nh - i]; } } return best; }
    case 14: return evenr(log((double)p));
    case 15: return evenr(log((double)p) / 2);
    case 16: return tab30[p % 30];
    case 17: return tab210[p % 210];
    case 18: return tab2310[p % 2310];
    case 19: { int k = (int)(p % 210) * 128 + (g1 / 2 < 127 ? g1 / 2 : 127); return tab210g[k] ? tab210g[k] : tab210[p % 210]; }
  }
  return 2;
}
static u64 predict(int r, u64 p, const int *h, int nh, int y) {
  int wk; u64 c = p + base(r, p, h, nh, &wk);
  if (y > 0) { int guard = 0; while (hasf(c, y) && guard++ < 100000) c += wk ? wheelgap(c, wk) : (c & 1 ? 2 : 1); }
  return c;
}
static void learn(int *tab, int M, const u64 *pr, int n, int withg) {
  // count gap frequencies per key, keep argmax (ties: smaller gap)
  int K = withg ? M * 128 : M, G = 512;
  int *cnt = calloc((size_t)K * G, sizeof(int));
  for (int i = 1; i + 1 < n; i++) {
    int g = (int)(pr[i + 1] - pr[i]); if (g >= G) continue;
    int key = (int)(pr[i] % M); if (withg) { int g1 = (int)(pr[i] - pr[i - 1]); key = key * 128 + (g1 / 2 < 127 ? g1 / 2 : 127); }
    cnt[(size_t)key * G + g]++;
  }
  for (int k = 0; k < K; k++) { int best = 0, bc = 0; for (int g = 1; g < G; g++) if (cnt[(size_t)k * G + g] > bc) { bc = cnt[(size_t)k * G + g]; best = g; } tab[k] = withg ? best : (best ? best : 2); }
  free(cnt);
}

int main(int argc, char **argv) {
  u64 X = strtoull(argv[1], 0, 10), T = strtoull(argv[2], 0, 10); int N = atoi(argv[3]);
  for (int i = 2; i <= 1000; i++) { int ok = 1; for (int j = 2; j * j <= i; j++) if (i % j == 0) { ok = 0; break; } if (ok) SMALL[nsmall++] = i; }
  u64 lo = X - T, hi = X + (u64)(N * 3.0 * log((double)X)) + 10000;
  size_t L = hi - lo + 1; unsigned char *comp = calloc(L, 1);
  u64 sq = (u64)sqrtl((long double)hi) + 2; unsigned char *bs = calloc(sq + 1, 1);
  for (u64 i = 2; i <= sq; i++) if (!bs[i]) { for (u64 j = i * i; j <= sq; j += i) bs[j] = 1; u64 s = (lo + i - 1) / i * i; if (s < i * i) s = i * i; for (u64 j = s; j <= hi; j += i) comp[j - lo] = 1; }
  u64 *pr = malloc(sizeof(u64) * (L / 2 + 10)); int np = 0;
  for (u64 n = lo < 2 ? 2 : lo; n <= hi; n++) if (!comp[n - lo]) pr[np++] = n;
  int first = 0; while (pr[first] < X) first++;
  if (np - first < N + 1) { fprintf(stderr, "window too small\n"); return 1; }
  tab30 = malloc(30 * sizeof(int)); tab210 = malloc(210 * sizeof(int)); tab2310 = malloc(2310 * sizeof(int)); tab210g = malloc(210 * 128 * sizeof(int));
  learn(tab30, 30, pr, first, 0); learn(tab210, 210, pr, first, 0); learn(tab2310, 2310, pr, first, 0); learn(tab210g, 210, pr, first, 1);
  int *gaps = malloc(sizeof(int) * np); for (int i = 1; i < np; i++) gaps[i] = (int)(pr[i] - pr[i - 1]);
  int H = 64; int *hb = malloc(sizeof(int) * (H + 4000));
  fprintf(stderr, "X=%llu train primes=%d test=%d\n", (unsigned long long)X, first, N);
  for (int r = 0; r < NR; r++) for (int yi = 0; yi < 5; yi++) {
    int y = YS[yi]; long hits = 0, sh = 0, lg = 0;
    for (int i = first; i < first + N; i++) {
      u64 c = predict(r, pr[i], gaps + i - H + 1, H, y);   // history = true gaps up to p
      if (c == pr[i + 1]) hits++; else if (c < pr[i + 1]) sh++; else lg++;
    }
    // own-power runs from the first 1000 test primes: primes in a row (incl. start), cap 2000
    double rs = 0; int rmax = 0; const int RUNS = 500, CAP = 500;
    for (int s = 0; s < RUNS; s++) {
      int i = first + s; memcpy(hb, gaps + i - H + 1, sizeof(int) * H); int nh = H;
      u64 p = pr[i]; int len = 1;
      while (len < CAP) { u64 c = predict(r, p, hb + nh - H, H, y); if (!isprime(c)) break; hb[nh++] = (int)(c - p); p = c; len++; if (nh >= H + 3900) { memmove(hb, hb + nh - H, sizeof(int) * H); nh = H; } }
      rs += len; if (len > rmax) rmax = len;
    }
    printf("%d,%d,%ld,%ld,%ld,%.3f,%d\n", r, y, hits, sh, lg, rs / RUNS, rmax);
    fflush(stdout);
  }
  return 0;
}
