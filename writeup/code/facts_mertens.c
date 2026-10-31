// facts_mertens.c: Mertens M(x) and Liouville L(x) for all x <= X, for FACTS.md section 2. Build: gcc -O3 -o facts_mertens facts_mertens.c -lm
// Run: ./facts_mertens 10000000000 (about 3 minutes). Segmented: each n is factored by the primes up to sqrt(X), and the leftover cofactor is one prime.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64; typedef int64_t i64;
int main(int argc, char **argv) {
  u64 X = strtoull(argv[1], 0, 10), R = (u64)sqrt((double)X) + 2;
  unsigned char *bs = calloc(R + 1, 1); u64 *bp = malloc(sizeof(u64) * R); long nb = 0;
  for (u64 i = 2; i <= R; i++) if (!bs[i]) { bp[nb++] = i; for (u64 j = i * i; j <= R; j += i) bs[j] = 1; }
  const u64 S = 1 << 20; u64 *prod = malloc(8 * S); signed char *mu = malloc(S); unsigned char *par = malloc(S);
  u64 watch[] = {906150257ULL, 906180359ULL, 906316571ULL, 906488079ULL, 906488080ULL, 7766842813ULL, 0};
  i64 M = 0, L = 0, Lmax = -1000000000; u64 Lmax_at = 0, Lfirst = 0, Llast = 0, Lposcount = 0, nextck = 10;
  double best = 0, bestneg = 0; u64 best_at = 0, bestneg_at = 0;
  for (u64 lo = 1; lo <= X; lo += S) {
    u64 hi = lo + S; if (hi > X + 1) hi = X + 1; u64 n = hi - lo;
    for (u64 i = 0; i < n; i++) prod[i] = 1; memset(mu, 1, n); memset(par, 0, n);
    for (long t = 0; t < nb; t++) { u64 p = bp[t]; if (p * p >= hi) break;
      u64 st = ((lo + p - 1) / p) * p;
      for (u64 j = st; j < hi; j += p) { prod[j - lo] *= p; mu[j - lo] = -mu[j - lo]; par[j - lo] ^= 1; }
      u64 pk = p * p, st2;
      if (pk < hi) { st2 = ((lo + pk - 1) / pk) * pk; for (u64 j = st2; j < hi; j += pk) mu[j - lo] = 0; }
      while (pk < hi) { st2 = ((lo + pk - 1) / pk) * pk;
        for (u64 j = st2; j < hi; j += pk) { prod[j - lo] *= p; par[j - lo] ^= 1; }
        if (pk > hi / p) break; pk *= p; }
    }
    for (u64 i = 0; i < n; i++) { u64 x = lo + i;
      int m = mu[i], pa = par[i];
      if (prod[i] < x) { m = -m; pa ^= 1; }
      M += m; L += pa ? -1 : 1;
      if (x >= 2 && L > 0) { Lposcount++; if (!Lfirst) Lfirst = x; Llast = x; }
      if (x >= 2 && L > Lmax) { Lmax = L; Lmax_at = x; }
      if (x >= 10000) { double r = (double)M * fabs((double)M) / (double)x;
        if (r > best) { best = r; best_at = x; } if (r < bestneg) { bestneg = r; bestneg_at = x; } }
      if (x == nextck) { printf("CK %llu M %lld L %lld\n", (unsigned long long)x, (long long)M, (long long)L); fflush(stdout); nextck *= 10; }
      for (int w = 0; watch[w]; w++) if (x == watch[w]) { printf("AT %llu M %lld L %lld\n", (unsigned long long)x, (long long)M, (long long)L); fflush(stdout); }
    }
  }
  printf("LFIRST %llu LLAST %llu LPOS %llu LMAX %lld at %llu (max over x >= 2)\n", (unsigned long long)Lfirst, (unsigned long long)Llast,
    (unsigned long long)Lposcount, (long long)Lmax, (unsigned long long)Lmax_at);
  printf("MMAX M/sqrt(x) = %.6f at %llu ; MMIN %.6f at %llu (x >= 10^4)\n", sqrt(best), (unsigned long long)best_at, -sqrt(-bestneg), (unsigned long long)bestneg_at);
}
