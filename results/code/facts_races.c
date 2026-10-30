// facts_races.c: Kronecker races and the mod 8 / mod 12 races, for FACTS.md section 1F. Build: gcc -O3 -o facts_races facts_races.c -lm
// Run: ./facts_races 140000000000 (about 8 minutes). For each k = +-1..+-64 it prints the first prime p with sum_{q<=p} (k/q) > 0.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64; typedef int64_t i64;
static int jacobi(i64 a, i64 n) { a %= n; if (a < 0) a += n; int r = 1;
  while (a) { while (!(a & 1)) { a >>= 1; i64 m = n & 7; if (m == 3 || m == 5) r = -r; }
    i64 t = a; a = n; n = t; if ((a & 3) == 3 && (n & 3) == 3) r = -r; a %= n; }
  return n == 1 ? r : 0; }
static int kron(i64 k, u64 p) { if (p == 2) { if (!(k & 1)) return 0; i64 m = ((k % 8) + 8) % 8; return (m == 1 || m == 7) ? 1 : -1; }
  return jacobi(k, (i64)p); }
#define NK 128
static i64 K[NK], S[NK]; static u64 F[NK]; static int per[NK]; static signed char *tab[NK]; static int active;
static i64 c8[8], c12[12]; static u64 f8[8], f12[12];
static void add(u64 p) {
  for (int i = 0; i < NK; i++) if (!F[i]) {
    int v = (p <= (u64)per[i] || p % 2 == 0 || (K[i] % (i64)p) == 0) ? kron(K[i], p) : tab[i][p % per[i]];
    S[i] += v; if (S[i] > 0) { F[i] = p; active--; }
  }
  c8[p % 8]++; c12[p % 12]++;
  for (int a = 3; a < 8; a += 2) if (!f8[a] && c8[1] > c8[a]) { f8[a] = p; printf("MOD8 1>%d first at %llu\n", a, (unsigned long long)p); fflush(stdout); }
  for (int a = 5; a < 12; a += 2) if (a != 9 && !f12[a] && c12[1] > c12[a]) { f12[a] = p; printf("MOD12 1>%d first at %llu\n", a, (unsigned long long)p); fflush(stdout); }
}
int main(int argc, char **argv) {
  u64 X = strtoull(argv[1], 0, 10);
  for (int i = 0; i < 64; i++) { K[i] = i + 1; K[64 + i] = -(i + 1); }
  active = NK;
  u64 R = (u64)sqrt((double)X) + 2; unsigned char *bs = calloc(R + 1, 1); u64 *bp = malloc(sizeof(u64) * R); long nb = 0;
  for (u64 i = 2; i <= R; i++) if (!bs[i]) { if (i > 2) bp[nb++] = i; for (u64 j = i * i; j <= R; j += i) bs[j] = 1; }
  for (int i = 0; i < NK; i++) { per[i] = 4 * (int)llabs(K[i]); tab[i] = calloc(per[i], 1);
    for (int r = 1; r < per[i]; r++) { u64 g = r, h = per[i]; while (h) { u64 t = g % h; g = h; h = t; } if (g != 1) continue;
      for (u64 q = r + per[i]; ; q += per[i]) { int pr = 1; for (u64 d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
        if (pr) { tab[i][r] = kron(K[i], q); break; } } } }
  add(2);
  const u64 SS = 1 << 22; unsigned char *seg = malloc(SS);
  for (u64 base = 3; base <= X; base += 2 * SS) {
    memset(seg, 0, SS); u64 top = base + 2 * SS;
    for (long t = 0; t < nb; t++) { u64 p = bp[t]; if (p * p >= top) break;
      u64 st = p * p; if (st < base) { st = ((base + p - 1) / p) * p; if (!(st & 1)) st += p; }
      for (u64 j = st; j < top; j += 2 * p) seg[(j - base) / 2] = 1; }
    for (u64 i = 0; i < SS; i++) { u64 p = base + 2 * i; if (p > X) break; if (!seg[i]) add(p); }
  }
  printf("POS"); for (int i = 0; i < 64; i++) printf(" %llu", (unsigned long long)F[i]); printf("\n");
  printf("NEG"); for (int i = 64; i < NK; i++) printf(" %llu", (unsigned long long)F[i]); printf("\n");
}
