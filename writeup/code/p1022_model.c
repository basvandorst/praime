/* Puzzle 1022: how many primes x(n) = (n-1)^n + n should there be?  (writeup/1022.md)

   Model: sieve every n exactly by each prime r <= R (x(n) mod r depends only on n mod r(r-1)),
   remove the algebraic class n = 2 (mod 3) (factor n^2 - n + 1) and every n with gcd(x(n), n+1) > 1
   or gcd(x(n), 2n-1) > 1 (the two divisor families of the write-up).  A survivor is prime with
   probability e^gamma ln R / ln x(n) (Mertens).  The same sieve and weight run on two sibling families
   whose prime indices are published, as a calibration:
     fam 0: x(n) = (n-1)^n + n
     fam 1: 3 n^n + 2        (OEIS A160360: 93, 605 in [30, 45000], next term > 45000)
     fam 2: (n-1) n^n + 1    (OEIS A108879: 93, 100, 1536, 2835 in [30, 2835])

   Build: gcc -O2 -o p1022_model p1022_model.c -lm
   Usage: ./p1022_model FAM LO HI R
   Prints the expected number of primes for LO <= n <= HI and, for fam 0, the cumulative expected
   count at checkpoints, the survivors by n mod 4 and the mean constant c = e^gamma ln R * [survivor]
   per decade.  ./p1022_model 0 30 13389 10000 takes 0.2 s; ./p1022_model 0 13390 10000000 10000
   about 3 minutes. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef unsigned __int128 u128;

static uint64_t powmod(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t r = 1 % m; b %= m;
    while (e) { if (e & 1) r = (u128)r * b % m; b = (u128)b * b % m; e >>= 1; }
    return r;
}
static uint64_t gcd(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }

/* f(n) mod r for the chosen family; the exponent is reduced to ((n-1) mod (r-1)) + 1 >= 1,
   which equals n mod (r-1) in the exponent and keeps 0^n = 0. */
static uint64_t fmod_r(int fam, uint64_t n, uint64_t r) {
    uint64_t e = ((n - 1) % (r - 1)) + 1;
    if (fam == 0) return (powmod(n - 1, e, r) + n) % r;
    if (fam == 1) return (3 * powmod(n, e, r) + 2) % r;
    return ((n - 1) % r * powmod(n, e, r) + 1) % r;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s FAM LO HI R\n", argv[0]); return 1; }
    int fam = atoi(argv[1]);
    long LO = atol(argv[2]), HI = atol(argv[3]), R = atol(argv[4]);
    char *comp = calloc(R + 1, 1);
    long *pr = malloc(sizeof(long) * (R + 1)); int np = 0;
    for (long i = 2; i <= R; i++) if (!comp[i]) { pr[np++] = i; for (long j = i * i; j <= R; j += i) comp[j] = 1; }
    double eg = exp(0.5772156649015329) * log((double)R), E = 0, wsum = 0;
    long cps[] = {13389, 20000, 30000, 50000, 100000, 200000, 300000, 500000, 700000,
                  1000000, 2000000, 5000000, 10000000, 100000000};
    int ci = 0; while (ci < 14 && cps[ci] < LO) ci++;
    long cat[5] = {0}, bymod4[4] = {0}, wn = 0, w0 = LO;   /* cat: 0 alg, 1 small r, 2 n+1, 3 2n-1, 4 survivor */
    for (long n = LO; n <= HI; n++) {
        int c = 4;
        if (fam == 0 && n % 3 == 2) c = 0;
        for (int k = 0; c == 4 && k < np; k++)
            if (fmod_r(fam, n, pr[k]) == 0) c = 1;
        if (fam == 0 && c == 4) {
            uint64_t m = n + 1;
            if (gcd((powmod(n - 1, n, m) + n) % m, m) > 1) c = 2;
        }
        if (fam == 0 && c == 4) {
            uint64_t m = 2 * n - 1;
            if (gcd((powmod(n - 1, n, m) + n) % m, m) > 1) c = 3;
        }
        cat[c]++;
        if (c == 4) {
            double lx = fam == 0 ? n * log((double)(n - 1))
                      : fam == 1 ? n * log((double)n) + log(3.0)
                                 : n * log((double)n) + log((double)(n - 1));
            E += eg / lx; wsum += eg; bymod4[n % 4]++;
        }
        wn++;
        if (fam == 0 && (n + 1 == w0 * 10 || n == HI)) {
            printf("  n in [%ld, %ld]: c = mean(e^gamma ln R * [survivor]) = %.4f\n", w0, n, wsum / wn);
            w0 = n + 1; wsum = 0; wn = 0;
        }
        if (fam == 0 && ci < 14 && n == cps[ci]) {
            printf("n <= %ld: cumulative expected primes from n = %ld: %.4f\n", n, LO, E);
            ci++;
        }
    }
    printf("fam %d, %ld <= n <= %ld, R = %ld: expected primes %.4f\n", fam, LO, HI, R, E);
    if (fam == 0)
        printf("  n = 2 mod 3: %ld, factor <= R: %ld, gcd with n+1: %ld, gcd with 2n-1: %ld, survivors: %ld "
               "(by n mod 4: %ld %ld %ld %ld)\n", cat[0], cat[1], cat[2], cat[3], cat[4],
               bymod4[0], bymod4[1], bymod4[2], bymod4[3]);
    else
        printf("  survivors: %ld\n", cat[4]);
    return 0;
}
