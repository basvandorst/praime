// Build: gcc -O2 -DVAR=2 -o var2 p377_search.c  (or -DVAR=3); run: ./var2 <k>  for an exhaustive search of k-digit n.
// Puzzle 377 variation #3: k-digit n, all digits 2..9, with n - 2*10^i prime for every i = 0..k-1 (base n composite).
// Exhaustive for one k: n = H*10^h + L. L buckets by L mod M (M = 3*7*13*17*19); for each H only the residues r
// with (H*10^h + r) - 2*10^i nonzero mod every q | M are visited. Then trial residues for q < 1000, then Miller-Rabin.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef VAR
#define VAR 3
#endif
#if VAR == 2
#define DLO 0
#define LA 1
#define LB 7
#define SG(n,x) ((n) + (x))
#else
#define DLO 2
#define LA 3
#define LB 9
#define SG(n,x) ((n) - (x))
#endif
static uint64_t lim0;
typedef unsigned __int128 u128;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (u128)a * b % m; }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) { uint64_t r = 1; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static int isprime(uint64_t n) {
  if (n < 2) return 0;
  static const uint64_t B[] = {2,3,5,7,11,13,17,19,23,29,31,37};
  for (int i = 0; i < 12; i++) { if (n % B[i] == 0) return n == B[i]; }
  uint64_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  for (int i = 0; i < 12; i++) {
    uint64_t x = powmod(B[i], d, n);
    if (x == 1 || x == n - 1) continue;
    int ok = 0; for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
    if (!ok) return 0;
  }
  return 1;  // deterministic for n < 3.3e24 with these bases
}
static int gen(int len, int lastrestrict, uint64_t *out) {
  lim0 = 1; for (int i = 1; i < len; i++) lim0 *= 10;   /* H must have a nonzero leading digit */   // all len-digit strings over 2..9 (last digit 3 or 9 if asked)
  int cnt = 0; uint64_t lim = 1; for (int i = 0; i < len; i++) lim *= 8;
  for (uint64_t c = 0; c < lim; c++) {
    uint64_t v = 0, t = c; int last = 0;
    for (int i = 0; i < len; i++) { int d = DLO + t % 8; t /= 8; if (i == len - 1) last = d; v = v * 10 + d; }
    if (lastrestrict && last != LA && last != LB) continue;
    if (!lastrestrict && len > 0 && v < lim0) continue;
    out[cnt++] = v;
  }
  return cnt;
}
int main(int argc, char **argv) {
  int k = atoi(argv[1]); int h = k / 2; if (h > 6) h = 6; int hl = k - h;   // L has h digits, H has hl digits
  uint64_t p10[20]; p10[0] = 1; for (int i = 1; i < 20; i++) p10[i] = p10[i-1] * 10;
  uint64_t nH = 1, nL = 1; for (int i = 0; i < hl; i++) nH *= 8; for (int i = 0; i < h; i++) nL *= 8;
  uint64_t *H = malloc(nH * 8), *L = malloc(nL * 8);
  nH = gen(hl, 0, H); nL = gen(h, 1, L);
  const int MQ[] = {3, 7, 13, 17, 19}; const int M = 3*7*13*17*19;
  // bad[r] = 1 if r (as n mod q for some q|M) hits n == 2*10^i mod q
  char *badM = calloc(M, 1);
  for (int r = 0; r < M; r++) for (int j = 0; j < 5; j++) for (int i = 0; i < k; i++)
    if ((uint64_t)r % MQ[j] == ((VAR == 2) ? (MQ[j] - 2 * p10[i] % MQ[j]) % MQ[j] : 2 * p10[i] % MQ[j])) badM[r] = 1;
  int *allowed = malloc(M * sizeof(int)); int na = 0; for (int r = 0; r < M; r++) if (!badM[r]) allowed[na++] = r;
  // bucket L by L mod M
  int *head = malloc(M * sizeof(int)), *next = malloc(nL * sizeof(int)); for (int r = 0; r < M; r++) head[r] = -1;
  for (uint64_t j = 0; j < nL; j++) { int r = L[j] % M; next[j] = head[r]; head[r] = j; }
  int primes[200], np = 0; for (int q = 23; q < 1000 && np < 200; q++) { int ok = 1; for (int d = 2; d * d <= q; d++) if (q % d == 0) { ok = 0; break; } if (ok) primes[np++] = q; }
  long long cand = 0, sol = 0;
  for (uint64_t a = 0; a < nH; a++) {
    uint64_t base = H[a] * p10[h]; int bM = base % M;
    for (int t = 0; t < na; t++) {
      int need = (allowed[t] - bM + M) % M;
      for (int j = head[need]; j >= 0; j = next[j]) {
        uint64_t n = base + L[j]; cand++;
        int bad = 0;
        for (int s = 0; s < np && !bad; s++) { int q = primes[s]; uint64_t nq = n % q; for (int i = 0; i < k; i++) { uint64_t t = (VAR == 2) ? (q - 2 * p10[i] % q) % q : 2 * p10[i] % q; if (nq == t && SG(n, 2 * p10[i]) != (uint64_t)q) { bad = 1; break; } } }
        if (bad) continue;
        int ok = 1; for (int i = 0; i < k && ok; i++) ok = isprime(SG(n, 2 * p10[i]));
        if (ok) { sol++; printf("k=%d solution %llu%s\n", k, (unsigned long long)n, isprime(n) ? " (prime base)" : ""); fflush(stdout); }
      }
    }
  }
  fprintf(stderr, "k=%d: %lld candidates after mod-M join, %lld solutions\n", k, cand, sol);
  return 0;
}
