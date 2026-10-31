// facts_gapfirst.c: first occurrence of every prime gap below X, for FACTS.md fact 3019. Build: gcc -O3 -o facts_gapfirst facts_gapfirst.c -lm
// Run: ./facts_gapfirst 140000000000 (about 8 minutes). Prints "g p" for each even gap g, p the prime before its first occurrence.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64;
#define MAXG 4096
int main(int argc, char **argv) {
  u64 X = strtoull(argv[1], 0, 10), R = (u64)sqrt((double)X) + 2, first[MAXG] = {0}, prev = 3;
  unsigned char *bs = calloc(R + 1, 1); u64 *bp = malloc(sizeof(u64) * R); long nb = 0;
  for (u64 i = 2; i <= R; i++) if (!bs[i]) { if (i > 2) bp[nb++] = i; for (u64 j = i * i; j <= R; j += i) bs[j] = 1; }
  const u64 S = 1 << 22; unsigned char *seg = malloc(S);
  for (u64 base = 5; base <= X; base += 2 * S) {
    memset(seg, 0, S); u64 top = base + 2 * S;
    for (long t = 0; t < nb; t++) { u64 p = bp[t]; if (p * p >= top) break;
      u64 st = p * p; if (st < base) { st = ((base + p - 1) / p) * p; if (!(st & 1)) st += p; }
      for (u64 j = st; j < top; j += 2 * p) seg[(j - base) / 2] = 1; }
    for (u64 i = 0; i < S; i++) { u64 p = base + 2 * i; if (p > X) break;
      if (!seg[i]) { u64 g = p - prev; if (g < MAXG && !first[g]) first[g] = prev; prev = p; } }
  }
  for (int g = 2; g < MAXG; g += 2) if (first[g]) printf("%d %llu\n", g, (unsigned long long)first[g]);
}
