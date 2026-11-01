/* Puzzle 207 (writeup/207.md): the digit-count map Phi of the inventory function on count vectors c[0..9] with sum <= 21 (every inventory sequence ends up with at most 21 digits).
   Phi(c)[e] = #{d : c_d > 0, d = e} + sum over d with c_d > 0 of (number of digits e in decimal(c_d)).
   Enumerates all vectors (stars and bars ranking), follows Phi, prints every cycle. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define D 10
#define S 21
static long long C[64][64];
static long long rank_(const int *c){          /* rank among vectors with sum <= S: use slack as 11th part, total = S */
  long long r = 0; int rem = S;
  for(int i=0;i<D;i++){ for(int x=0;x<c[i];x++) r += C[rem - x + (D - i) - 1][(D - i) - 1]; rem -= c[i]; }
  return r;
}
static void unrank(long long r, int *c){
  int rem = S;
  for(int i=0;i<D;i++){ int x=0; while(1){ long long cnt = C[rem - x + (D - i) - 1][(D - i) - 1]; if(r < cnt) break; r -= cnt; x++; } c[i]=x; rem -= x; }
}
static void phi(const int *c, int *o){
  memset(o,0,sizeof(int)*D);
  for(int d=0;d<D;d++) if(c[d]>0){ o[d]++; int k=c[d]; if(k>=10){ o[k/10]++; o[k%10]++; } else o[k]++; }
}
int main(){
  for(int n=0;n<64;n++){ C[n][0]=1; for(int k=1;k<=n;k++) C[n][k]=C[n-1][k-1]+(k<=n-1?C[n-1][k]:0); }
  long long N = C[S + D][D];                   /* vectors with sum <= S */
  fprintf(stderr, "states %lld\n", N);
  unsigned char *col = calloc(N, 1);           /* 0 new, 1 on stack, 2 done */
  int c[D], o[D]; long long *path = malloc(sizeof(long long) * 4096);
  int maxper = 0, ncyc = 0;
  for(long long s=0;s<N;s++){
    if(col[s]) continue;
    int len = 0; long long x = s;
    while(1){
      if(col[x] == 2) break;
      if(col[x] == 1){                           /* found a cycle: from x along path */
        int i = len - 1; while(path[i] != x) i--;
        int per = len - i; ncyc++; if(per > maxper) maxper = per;
        printf("cycle period %d:", per);
        for(int j=i;j<len;j++){ unrank(path[j], c); printf(" ["); for(int d=0;d<D;d++) printf("%d%s", c[d], d<D-1?",":""); printf("]"); }
        printf("\n");
        break;
      }
      col[x] = 1; path[len++] = x;
      unrank(x, c); phi(c, o);
      int sum=0; for(int d=0;d<D;d++) sum += o[d];
      if(sum > S){ break; }                      /* leaves the space (not part of any cycle) */
      x = rank_(o);
    }
    for(int j=0;j<len;j++) col[path[j]] = 2;
  }
  fprintf(stderr, "cycles %d, max count-period %d\n", ncyc, maxper);
  return 0;
}
