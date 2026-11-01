/* Puzzle 207 (writeup/207.md): for each cycle of the count map (cycles.txt), the order of first appearance of the digits evolves by a
   deterministic map. For every permutation sigma of the support of the first count vector, follow one full turn of
   the count cycle (P' steps) and get F(sigma). The period of a number cycle is P' * (cycle length of F).
   Prints, per count cycle, the cycle-length distribution of F, and the maximal number period overall. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int fact[11];
static int lehmer(const int *perm, int k){ int r=0; for(int i=0;i<k;i++){ int s=0; for(int j=i+1;j<k;j++) if(perm[j]<perm[i]) s++; r = r*(k-i) + s; } return r; }
static void unlehmer(int r, int k, int *perm){ int code[10]; for(int i=k-1;i>=0;i--){ code[i] = r % (k-i); r /= (k-i); }
  int used[10]={0}; for(int i=0;i<k;i++){ int c=code[i]; for(int j=0;j<k;j++) if(!used[j]){ if(c==0){ perm[i]=j; used[j]=1; break; } c--; } } }
int main(){
  fact[0]=1; for(int i=1;i<=10;i++) fact[i]=fact[i-1]*i;
  FILE *f = fopen("cycles.txt","r"); char line[8192]; int maxP = 0, ncyc = 0;
  while(fgets(line, sizeof line, f)){
    int per; sscanf(line, "cycle period %d:", &per);
    int c[3][10]; char *p = strchr(line, '['); for(int i=0;i<per;i++){ for(int d=0;d<10;d++){ c[i][d] = strtol(p+1, &p, 10); } p = strchr(p, '['); }
    ncyc++;
    int sup[10], k=0; for(int d=0;d<10;d++) if(c[0][d]>0) sup[k++]=d;      /* sorted support of c0 */
    if(k == 0){ printf("cycle %d (P'=%d): empty\n", ncyc, per); continue; }
    int n = fact[k]; int *F = malloc(sizeof(int)*n);
    for(int r=0;r<n;r++){
      int idx[10]; unlehmer(r, k, idx);
      int sigma[10], sl = k; for(int i=0;i<k;i++) sigma[i] = sup[idx[i]];
      for(int step=0; step<per; step++){
        int *cc = c[step]; char s[64]; int len = 0;
        for(int i=0;i<sl;i++){ int d = sigma[i], cnt = cc[d]; if(cnt>=10){ s[len++]='0'+cnt/10; s[len++]='0'+cnt%10; } else s[len++]='0'+cnt; s[len++]='0'+d; }
        int seen[10]={0}, nc[10]={0}; sl = 0;
        for(int i=0;i<len;i++){ int d = s[i]-'0'; nc[d]++; if(!seen[d]){ seen[d]=1; sigma[sl++]=d; } }
        int *cn = c[(step+1)%per]; for(int d=0;d<10;d++) if(nc[d] != cn[d]){ printf("count mismatch\n"); return 1; }
      }
      /* sigma is now a permutation of the support of c0 */
      int pos[10]; for(int i=0;i<k;i++) pos[sup[i]] = i;
      int perm[10]; for(int i=0;i<sl;i++) perm[i] = pos[sigma[i]];
      if(sl != k){ printf("support mismatch\n"); return 1; }
      F[r] = lehmer(perm, k);
    }
    /* cycle lengths of F */
    char *col = calloc(n,1); int *stk = malloc(sizeof(int)*n); int maxL = 0; long hist[64]={0};
    for(int s0=0;s0<n;s0++){ if(col[s0]) continue; int len=0, x=s0;
      while(!col[x]){ col[x]=1; stk[len++]=x; x=F[x]; }
      if(col[x]==1){ int i=len-1; while(stk[i]!=x) i--; int L = len-i; if(L>maxL) maxL=L; if(L<64) hist[L]++; }
      for(int j=0;j<len;j++) col[stk[j]]=2; }
    int P = per*maxL; if(P>maxP) maxP=P;
    printf("count cycle %3d: P'=%d support %2d: F-cycle lengths", ncyc, per, k);
    for(int L=1;L<64;L++) if(hist[L]) printf(" %d(x%ld)", L, hist[L]);
    printf("  -> max number period %d\n", P);
    free(F); free(col); free(stk);
  }
  printf("MAX number-cycle period over all count cycles: %d\n", maxP);
  return 0;
}
