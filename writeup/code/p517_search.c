/* Puzzle 517: all x <= N with phi(x) + sigma(x+1) = 3x.  phi and sigma by sieves (32-bit).
   Usage: p517_search N  (N <= 2e8) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv) {
    uint32_t N = (uint32_t)atoll(argv[1]);
    uint32_t *phi = malloc((size_t)(N + 2) * 4), *sig = calloc(N + 2, 4);
    for (uint32_t i = 0; i <= N + 1; i++) phi[i] = i;
    for (uint32_t p = 2; p <= N + 1; p++) if (phi[p] == p) for (uint32_t m = p; m <= N + 1; m += p) phi[m] -= phi[m] / p;
    for (uint32_t d = 1; d <= N + 1; d++) for (uint32_t m = d; m <= N + 1; m += d) sig[m] += d;
    long cnt = 0;
    for (uint32_t x = 1; x <= N; x++) if ((uint64_t)phi[x] + sig[x + 1] == 3ULL * x) { cnt++; printf("x = %u  (x+1 = %u)\n", x, x + 1); }
    printf("N = %u: %ld solutions\n", N, cnt);
    return 0;
}
