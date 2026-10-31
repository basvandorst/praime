// Conjecture 21: least N for each R(N). Build: gcc -O3 -o conj21 conj-21.c   Run: ./conj21 100000000 12  (about 1 minute)
// Conjecture 21: R(N) = least r such that N is a signed sum of r consecutive primes, all <= N.
// Marks best[N] for N <= M using blocks of r <= RMAX primes; sign patterns enumerated by Gray code.
// Prints the least N for each R and how many N <= M have each R.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv){
  int64_t M = atoll(argv[1]); int RMAX = atoi(argv[2]);
  char *comp = calloc(M + 1, 1); uint32_t *P = malloc(sizeof(uint32_t) * (M / 10 + 1000)); int64_t np = 0;
  for (int64_t i = 2; i <= M; i++) { if (!comp[i]) { P[np++] = i; for (int64_t j = i * i; j <= M; j += i) comp[j] = 1; } }
  free(comp);
  unsigned char *best = malloc(M + 1); for (int64_t i = 0; i <= M; i++) best[i] = 255;
  for (int r = 1; r <= RMAX; r++) {
    for (int64_t s = 0; s + r <= np; s++) {
      int64_t pmax = P[s + r - 1];
      if (pmax > M) break;
      // start with all signs +, then walk the Gray code (2^r patterns)
      int64_t sum = 0; for (int j = 0; j < r; j++) sum += P[s + j];
      if (sum >= pmax && sum <= M && best[sum] > r) best[sum] = r;
      unsigned int sign = 0;   // bit j set = element j negative
      for (uint32_t g = 1; g < (1u << r); g++) {
        int j = __builtin_ctz(g);           // bit that flips in the Gray code
        sign ^= 1u << j;
        if (sign >> j & 1) sum -= 2 * (int64_t)P[s + j]; else sum += 2 * (int64_t)P[s + j];
        if (sum >= pmax && sum <= M && best[sum] > r) best[sum] = r;
      }
    }
    fprintf(stderr, "r=%d done\n", r);
  }
  int64_t cnt[256] = {0}, first[256]; for (int i = 0; i < 256; i++) first[i] = -1;
  for (int64_t n = 2; n <= M; n++) { cnt[best[n]]++; if (first[best[n]] < 0) first[best[n]] = n; }
  for (int r = 1; r <= RMAX; r++) printf("R=%d least N=%lld count=%lld\n", r, (long long)first[r], (long long)cnt[r]);
  printf("unreached (R > %d): count=%lld, first=%lld", RMAX, (long long)cnt[255], (long long)first[255]);
  int shown = 0; for (int64_t n = 5; n <= M && shown < 10; n++) if (best[n] == 255) { printf(" %lld", (long long)n); shown++; }
  printf("\n"); return 0;
}
