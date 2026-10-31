// Count Euler-Rodriguez chains of exactly k consecutive primes: forward (gaps 2,4,...,2(k-1)) and
// mirrored (gaps 2(k-1),...,4,2), for FACTS.md fact 4010. Build: gcc -O3 -o facts_chains facts_chains.c -lm
// Run: ./facts_chains 10000000000   (about 13 s)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64;
#define GB 64
static u64 fwd[GB], rev[GB], ffirst[GB], rfirst[GB], gaps[GB], ng = 0, prevp = 0, pbuf[GB];
static int L = 0; static u64 A = 0;
static void add(u64 p) {
  if (prevp) {
    u64 g = p - prevp; gaps[ng % GB] = g; pbuf[ng % GB] = prevp; ng++;
    if (g == (u64)(2 * L)) L++;
    else { if (L >= 2) { fwd[L]++; if (!ffirst[L]) ffirst[L] = A; }
           if (g == 2) { L = 2; A = prevp; } else { L = 1; A = p; } }
    if (g == 2) {                                   // mirrored chain ending with this twin gap
      int k = 2;
      while (k < GB - 1 && (u64)k <= ng && gaps[(ng - k) % GB] == (u64)(2 * k)) k++;
      rev[k]++; if (!rfirst[k]) rfirst[k] = pbuf[(ng - (k - 1)) % GB];
    }
  }
  prevp = p;
}
int main(int argc, char **argv) {
  u64 X = strtoull(argv[1], 0, 10);
  u64 R = (u64)sqrt((double)X) + 2; char *bs = calloc(R + 1, 1); u64 *bp = malloc(sizeof(u64) * R); long nb = 0;
  for (u64 i = 2; i <= R; i++) if (!bs[i]) { if (i > 2) bp[nb++] = i; for (u64 j = i * i; j <= R; j += i) bs[j] = 1; }
  const u64 S = 1 << 22; unsigned char *seg = malloc(S);
  for (u64 base = 3; base <= X; base += 2 * S) {
    memset(seg, 0, S); u64 top = base + 2 * S;
    for (long t = 0; t < nb; t++) { u64 p = bp[t]; if (p * p >= top) break;
      u64 st = ((base + p - 1) / p) * p; if (st < p * p) st = p * p; if (st % 2 == 0) st += p;
      for (u64 j = st; j < top; j += 2 * p) seg[(j - base) / 2] = 1; }
    for (u64 i = 0; i < S; i++) { u64 p = base + 2 * i; if (p > X) break; if (!seg[i]) add(p); }
  }
  for (int k = 2; k < 20; k++) if (fwd[k] || rev[k])
    printf("k=%2d forward %12llu first %14llu   mirrored %12llu first %14llu\n", k, (unsigned long long)fwd[k],
      (unsigned long long)ffirst[k], (unsigned long long)rev[k], (unsigned long long)rfirst[k]);
  return 0;
}
