// Puzzle 896, set A, second method. Build: gcc -O2 -o p896_scs_bfs p896_scs_bfs.c ; run: ./p896_scs_bfs 35 (18 s: not reachable) / 36
// Independent check: breadth-first layers of position vectors (all words advanced greedily by the chosen digit).
// Layer t = all states reachable with t digits, deduplicated, keeping only states with t + LB <= L, where
// LB = sum over digits of the largest remaining count of that digit (digit-count bound only).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define K 10
static const char *W[K]={"4242064951","5519385407","6886477529","3816749731","4700006089","9650560979","4371455527","5419320281","8761823531","3454277693"};
static int cnt[K][11][10];
static void dec(uint64_t k, int *p){ for(int i=K-1;i>=0;i--){ p[i]=k%11; k/=11; } }
static uint64_t enc(const int *p){ uint64_t k=0; for(int i=0;i<K;i++) k=k*11+p[i]; return k; }
static int lbc(const int *p){ int s=0; for(int d=0;d<10;d++){ int m=0; for(int i=0;i<K;i++) if(cnt[i][p[i]][d]>m) m=cnt[i][p[i]][d]; s+=m; } return s; }
static int cmp(const void*a,const void*b){ uint64_t x=*(uint64_t*)a,y=*(uint64_t*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){ int L=atoi(argv[1]);
  for(int i=0;i<K;i++) for(int a=10;a>=0;a--) for(int d=0;d<10;d++) cnt[i][a][d]= a==10?0: cnt[i][a+1][d] + (W[i][a]-'0'==d);
  size_t cap=1<<27; uint64_t *cur=malloc(8*cap), *nxt=malloc(8*cap); size_t nc=1, nn; cur[0]=0;
  uint64_t goal; { int p[K]; for(int i=0;i<K;i++) p[i]=10; goal=enc(p); }
  for(int t=0;t<=L;t++){
    for(size_t s=0;s<nc;s++) if(cur[s]==goal){ printf("goal reached with %d digits\n",t); return 0; }
    if(t==L) break;
    nn=0;
    for(size_t s=0;s<nc;s++){ int p[K]; dec(cur[s],p);
      int seen[10]={0}; for(int i=0;i<K;i++) if(p[i]<10) seen[W[i][p[i]]-'0']=1;
      for(int d=0;d<10;d++) if(seen[d]){ int q[K]; for(int i=0;i<K;i++) q[i]=(p[i]<10&&W[i][p[i]]-'0'==d)?p[i]+1:p[i];
        if(t+1+lbc(q)<=L){ if(nn==cap){ fprintf(stderr,"overflow\n"); return 1; } nxt[nn++]=enc(q); } } }
    qsort(nxt,nn,8,cmp); size_t u=0; for(size_t s=0;s<nn;s++) if(u==0||nxt[s]!=nxt[u-1]) nxt[u++]=nxt[s];
    uint64_t *tmp=cur; cur=nxt; nxt=tmp; nc=u;
    fprintf(stderr,"layer %d: %zu states\n",t+1,nc);
  }
  printf("goal not reachable with %d digits\n",L); return 0; }
