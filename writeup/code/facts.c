// facts.c: statistics of the primes below X, for FACTS.md. Build: gcc -O3 -o facts facts.c -lm
// Run: ./facts 140000000000   (about 8 minutes; output lines are tagged, see the printf calls)
// One pass over the primes in increasing order (segmented sieve of Eratosthenes over odd numbers).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64; typedef int64_t i64;

#define MAXG 2000
#define B 128                                   // ring buffer of the last B primes
static u64 X, buf[B], cnt = 0, prevp = 0, twins = 0, ck[16]; static int nck = 0, ci = 0;
static u64 hist[MAXG], m3[3], m4[4], m10[10], maxg = 0;
static double best_merit = 0, best_cramer = 0, best_sqrt = 0;
static int  runR[5];                            // run length on a polynomial of degree d (d = 1..4)
static u64 poly_first[5][B];                     // first start of a run of length R on degree d
static int  erL = 0; static u64 erA = 0, er_first[64];  // Euler-Rodriguez chains: gaps 2, 4, 6, ...
static int  shiu_len[3]; static u64 shiu_first[3][10][B]; static const int MODS[3] = {3, 4, 10};
static u64 minspan[20], span_first[20];
static u64 lo_trans[10][10];                     // last-digit transitions among the first 10^8 primes
static i64 d4 = 0, d3 = 0, d3min = 0; static u64 d3min_at = 0, d4_regions = 0, d4_neg = 0, d3_first_neg = 0;

static inline u64 back(int k) { return buf[(cnt - 1 - k) % B]; }   // k = 0: the newest prime

static void checkpoint(u64 x) {
  printf("PI %llu %llu twins %llu mod4 %llu %llu mod3 %llu %llu mod10 %llu %llu %llu %llu d4 %lld\n",
    (unsigned long long)x, (unsigned long long)cnt, (unsigned long long)twins,
    (unsigned long long)m4[1], (unsigned long long)m4[3], (unsigned long long)m3[1], (unsigned long long)m3[2],
    (unsigned long long)m10[1], (unsigned long long)m10[3], (unsigned long long)m10[7], (unsigned long long)m10[9],
    (long long)d4);
  if (x >= 1000000) {                            // the most common gaps so far
    int top[6] = {0}; for (int t = 0; t < 6; t++) { u64 b = 0; for (int g = 2; g < MAXG; g += 2) {
      int used = 0; for (int s = 0; s < t; s++) if (top[s] == g) used = 1;
      if (!used && hist[g] > b) { b = hist[g]; top[t] = g; } } }
    printf("TOPGAPS %llu", (unsigned long long)x);
    for (int t = 0; t < 6; t++) printf(" %d:%llu", top[t], (unsigned long long)hist[top[t]]);
    printf("\n");
  }
}

static void add(u64 p) {
  while (ci < nck && p > ck[ci]) checkpoint(ck[ci++]);
  buf[cnt % B] = p; cnt++;
  m3[p % 3]++; m4[p % 4]++; m10[p % 10]++;
  if (p > 2) { if (p % 4 == 1) d4--; else d4++; }
  if (p > 3) { if (p % 3 == 1) d3--; else d3++; }
  if (p % 4 == 1 && d4 == -1) { d4_regions++; if (d4_regions <= 12) printf("D4NEG %llu\n", (unsigned long long)p); }
  if (d4 < 0) d4_neg++;
  if (d3 < 0 && !d3_first_neg) d3_first_neg = p;
  if (p > 100000 && (d3 < d3min || d3min_at == 0)) { d3min = d3; d3min_at = p; }
  if (prevp) {
    u64 g = p - prevp;
    if (g == 2) twins++;
    if (g < MAXG) hist[g]++;
    if (g > maxg) { maxg = g; printf("MAXGAP %llu after %llu\n", (unsigned long long)g, (unsigned long long)prevp); }
    if (prevp >= 11) {
      double L = log((double)prevp), m = g / L, c = g / (L * L), s = g / sqrt((double)prevp);
      if (m > best_merit) { best_merit = m; printf("MERIT %.4f gap %llu after %llu\n", m, (unsigned long long)g, (unsigned long long)prevp); }
      if (c > best_cramer) { best_cramer = c; printf("CRAMER %.4f gap %llu after %llu\n", c, (unsigned long long)g, (unsigned long long)prevp); }
      if (s > best_sqrt) { best_sqrt = s; printf("SQRT %.4f gap %llu after %llu\n", s, (unsigned long long)g, (unsigned long long)prevp); }
    }
    // Euler-Rodriguez chain: A, A+2, A+6, A+12, ... (gaps 2, 4, 6, ...); record the first chain of exactly L primes
    if (g == (u64)(2 * erL)) erL++;
    else {
      if (erL >= 2 && erL < 64 && !er_first[erL]) { er_first[erL] = erA; printf("ER %d %llu\n", erL, (unsigned long long)erA); }
      if (g == 2) { erL = 2; erA = prevp; } else { erL = 1; erA = p; }
    }
    if (prevp >= 7 && cnt <= 100000000) lo_trans[prevp % 10][p % 10]++;
  } else { erL = 1; erA = p; }
  // runs of consecutive primes on a polynomial of degree d: (d+1)-th differences vanish
  static const i64 C[6][7] = {{1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},{1,5,10,10,5,1}};
  for (int d = 1; d <= 4; d++) {
    if (cnt < (u64)(d + 2) || back(d + 1) < 3) { runR[d] = (int)(cnt < (u64)(d + 1) ? cnt : (u64)(d + 1)); continue; }
    i64 s = 0; for (int j = 0; j <= d + 1; j++) s += ((j & 1) ? -1 : 1) * C[d + 1][j] * (i64)back(j);
    if (s == 0) runR[d]++; else runR[d] = d + 1;
    int R = runR[d];
    if (R >= d + 2 && R < B && !poly_first[d][R]) { poly_first[d][R] = back(R - 1);
      printf("POLY d=%d run %d from %llu\n", d, R, (unsigned long long)back(R - 1)); }
  }
  // strings of consecutive primes in one residue class (Shiu)
  for (int t = 0; t < 3; t++) {
    int m = MODS[t];
    if (cnt >= 2 && p > 5 && back(1) > 5 && back(1) % m == p % m) shiu_len[t]++; else shiu_len[t] = 1;
    int L = shiu_len[t];
    if (L >= 2 && L < B && !shiu_first[t][p % m][L]) { shiu_first[t][p % m][L] = back(L - 1);
      if (L >= 8) printf("SHIU mod %d class %llu run %d from %llu\n", m, (unsigned long long)(p % m), L, (unsigned long long)back(L - 1)); }
  }
  // tightest clusters: least span of k consecutive primes, first prime > k
  for (int k = 2; k < 20 && (u64)k <= cnt; k++) {
    u64 a = back(k - 1); if (a <= (u64)k) continue;
    u64 sp = p - a; if (!minspan[k] || sp < minspan[k]) { minspan[k] = sp; span_first[k] = a; }
  }
  prevp = p;
}

int main(int argc, char **argv) {
  X = strtoull(argv[1], 0, 10);
  for (u64 x = 10; x <= X; x *= 10) ck[nck++] = x;
  if (ck[nck - 1] != X) ck[nck++] = X;
  u64 R = (u64)sqrt((double)X) + 2; char *bs = calloc(R + 1, 1); u64 *bp = malloc(sizeof(u64) * R); long nb = 0;
  for (u64 i = 2; i <= R; i++) if (!bs[i]) { if (i > 2) bp[nb++] = i; for (u64 j = i * i; j <= R; j += i) bs[j] = 1; }
  const u64 S = 1 << 22; unsigned char *seg = malloc(S);
  add(2);
  for (u64 base = 3; base <= X; base += 2 * S) {
    memset(seg, 0, S); u64 top = base + 2 * S;
    for (long t = 0; t < nb; t++) { u64 p = bp[t]; if (p * p >= top) break;
      u64 st = ((base + p - 1) / p) * p; if (st < p * p) st = p * p; if (st % 2 == 0) st += p;
      for (u64 j = st; j < top; j += 2 * p) seg[(j - base) / 2] = 1; }
    for (u64 i = 0; i < S; i++) { u64 p = base + 2 * i; if (p > X) break; if (!seg[i]) add(p); }
  }
  while (ci < nck) checkpoint(ck[ci++]);
  printf("END primes %llu last %llu maxgap %llu\n", (unsigned long long)cnt, (unsigned long long)prevp, (unsigned long long)maxg);
  printf("D4 regions where pi(x;4,1) > pi(x;4,3): %llu; primes inside them: %llu\n", (unsigned long long)d4_regions, (unsigned long long)d4_neg);
  printf("D3 first negative: %llu; min of pi(x;3,2)-pi(x;3,1) for p > 1e5: %lld at %llu\n", (unsigned long long)d3_first_neg, (long long)d3min, (unsigned long long)d3min_at);
  for (int k = 2; k < 20; k++) printf("SPAN k=%d min %llu first %llu\n", k, (unsigned long long)minspan[k], (unsigned long long)span_first[k]);
  for (int t = 0; t < 3; t++) for (int c = 0; c < MODS[t]; c++) { int L = 0; for (int l = 2; l < B; l++) if (shiu_first[t][c][l]) L = l;
      if (L) printf("SHIUMAX mod %d class %d longest %d from %llu\n", MODS[t], c, L, (unsigned long long)shiu_first[t][c][L]); }
  printf("LO first 1e8 primes (from 7), rows = last digit of p_n, cols = last digit of p_n+1 (1 3 7 9):\n");
  int D[4] = {1, 3, 7, 9};
  for (int a = 0; a < 4; a++) { printf("LO %d:", D[a]); for (int b = 0; b < 4; b++) printf(" %llu", (unsigned long long)lo_trans[D[a]][D[b]]); printf("\n"); }
  printf("HIST"); for (int g = 1; g < MAXG; g++) if (hist[g]) printf(" %d:%llu", g, (unsigned long long)hist[g]); printf("\n");
  return 0;
}
