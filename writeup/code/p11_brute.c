// Independent check (different method): for every start prime, count how many consecutive gaps stay pairwise distinct
// (fresh bitset per start). Prints the first start reaching each N. Plain Eratosthenes per segment, all integers.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64;
int main(int argc, char **argv) {
  u64 LO0 = strtoull(argv[1], 0, 10), XMAX = strtoull(argv[2], 0, 10);
  const u64 SEG = 50000000ull, PAD = 200000;
  u64 r = (u64)sqrtl((long double)(XMAX + PAD)) + 2; char *bs = calloc(r + 1, 1); u64 *bp = malloc(sizeof(u64) * r / 2); int nb = 0;
  for (u64 i = 2; i <= r; i++) if (!bs[i]) { bp[nb++] = i; for (u64 j = i * i; j <= r; j += i) bs[j] = 1; }
  char *s = malloc(SEG + PAD); u64 *P = malloc(sizeof(u64) * (SEG + PAD) / 2);
  static u64 first[64];
  for (u64 lo = LO0 < 2 ? 2 : LO0; lo < XMAX; lo += SEG) {
    u64 L = SEG + PAD; memset(s, 1, L);
    for (int k = 0; k < nb; k++) { u64 p = bp[k]; if (p * p >= lo + L) break;
      u64 st = (lo + p - 1) / p * p; if (st < p * p) st = p * p; for (u64 j = st; j < lo + L; j += p) s[j - lo] = 0; }
    long n = 0; for (u64 i = 0; i < L; i++) if (s[i] && lo + i >= 3) P[n++] = lo + i;
    for (long a = 0; a + 1 < n && P[a] < lo + SEG && P[a] < XMAX; a++) {
      static long stamp[2048]; static long gen = 0; gen++; int len = 0;
      for (long b = a; b + 1 < n && len < 63; b++) { u64 g = P[b + 1] - P[b]; if (g >= 2048 || stamp[g] == gen) break; stamp[g] = gen; len++; }
      for (int N = 2; N <= len; N++) if (!first[N]) first[N] = P[a];
    }
  }
  for (int N = 2; N < 64; N++) if (first[N]) printf("D %d %llu\n", N, (unsigned long long)first[N]);
  return 0;
}
