/* mersieve build timing: JS-style loops in C.
   For every prime p <= 500000 find the smallest k <= K with 2^p == 1 mod q, q = 2kp+1, q < 2^p-1.
   Modes (all numbers are doubles, as in JavaScript, except mode 3):
     0  r = (a*b) % q via fmod                      exact only while a*b < 2^53 (q < 9.49e7, K <= 94)
     1  r = a*b - q*floor(a*b/q), one correction    same range
     2  split mulmod: ((a*(b>>16) % q)*65536 + a*(b&65535)) % q, two fmods, exact for q < 2^30
     3  uint64 integer %, reference
     4  imul trick (JS: Math.imul): qh = floor(a*b/q) in doubles, r = (imul(a,b) - imul(qh,q)) | 0,
        one correction; exact for q < 2^30 (all k <= 1000 here)
   filter=1 skips k whose q is not +-1 mod 8 (they can never divide); filter=2 also skips q with a prime factor 3..13.
   Squaring uses the mulmod; "times 2" is r+=r; if(r>=q) r-=q.
   Usage: ms_speed K mode filter  -> prints time and writes speed_first_K_mode_filter.txt  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

static int mode;
static inline double mulmod(double a, double b, double q) {
    if (mode == 0) return fmod(a * b, q);
    if (mode == 1) {
        double t = a * b; double r = t - q * floor(t / q);
        if (r < 0) r += q; else if (r >= q) r -= q;
        return r;
    }
    if (mode == 2) {
        double bh = floor(b / 65536.0), bl = b - bh * 65536.0;
        double x = fmod(a * bh, q) * 65536.0 + a * bl;
        return fmod(x, q);
    }
    if (mode == 3) return (double)(((uint64_t)a * (uint64_t)b) % (uint64_t)q);
    {
        uint32_t ua = (uint32_t)a, ub = (uint32_t)b, uq = (uint32_t)q;
        uint32_t qh = (uint32_t)floor(a * b / q);
        int32_t r = (int32_t)(ua * ub - qh * uq);   /* low 32 bits, like Math.imul(...) - Math.imul(...) | 0 */
        if (r < 0) r += (int32_t)uq; else if (r >= (int32_t)uq) r -= (int32_t)uq;
        return (double)r;
    }
}
static int powtest(int p, double q) {   /* 2^p mod q == 1 ? */
    int nb = 31 - __builtin_clz((unsigned)p);
    double r = 1;
    for (int b = nb; b >= 0; b--) {
        r = mulmod(r, r, q);
        if ((p >> b) & 1) { r += r; if (r >= q) r -= q; }
    }
    return r == 1;
}
int main(int argc, char **argv) {
    int K = atoi(argv[1]); mode = atoi(argv[2]); int filter = atoi(argv[3]);
    const int N = 500000;
    char *s = calloc(N + 1, 1);
    for (int i = 2; i <= N; i++) s[i] = 1;
    for (int i = 2; (long)i * i <= N; i++) if (s[i]) for (int j = i * i; j <= N; j += i) s[j] = 0;
    int *res = calloc(N + 1, sizeof(int));
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    long tests = 0;
    for (int p = 2; p <= N; p++) {
        if (!s[p]) continue;
        int found = 0;
        for (int k = 1; k <= K; k++) {
            double q = 2.0 * k * p + 1;
            if (filter) { int m = (int)fmod(q, 8.0); if (m != 1 && m != 7) continue; }
            if (filter == 2) {   /* also skip q with a factor 3..13 (q != that prime): the smallest dividing q is prime */
                static const int sp[5] = {3, 5, 7, 11, 13}; int skip = 0;
                for (int j = 0; j < 5; j++) { double r = fmod(q, (double)sp[j]); if (r == 0 && q != sp[j]) { skip = 1; break; } }
                if (skip) continue;
            }
            if (p < 64 && q >= ldexp(1.0, p) - 1) break;   /* q must be a proper factor */
            tests++;
            if (powtest(p, q)) { found = k; break; }
        }
        res[p] = found;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double dt = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    printf("K=%d mode=%d filter=%d: %ld modpows, %.3f s, %.1f ns/modpow\n", K, mode, filter, tests, dt, 1e9 * dt / tests);
    char fn[64]; sprintf(fn, "speed_first_%d_%d_%d.txt", K, mode, filter);
    FILE *f = fopen(fn, "w");
    for (int p = 2; p <= N; p++) if (s[p]) fprintf(f, "%d %d\n", p, res[p]);
    fclose(f);
    return 0;
}
