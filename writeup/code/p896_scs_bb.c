// Puzzle 896, set A. Build: gcc -O2 -o p896_scs_bb p896_scs_bb.c ; run: ./p896_scs_bb 35  (21 s: NONE) / ./p896_scs_bb 36 (FOUND)
// Exact shortest common supersequence decision for K words of length <= 16: is there a string of length <= L?
// DFS over position vectors with lower bound max(pairwise SCS of suffixes, sum_d max_i count_d(suffix_i)),
// transposition table storing the smallest depth at which a state was expanded.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define K 10
static char w[K][17]; static int len[K];
static unsigned char pair[K][K][17][17];   // SCS length of suffixes w_i[a..], w_j[b..]
static unsigned char cnt[K][17][10];        // count of digit d in suffix w_i[a..]
static int L; static long long nodes=0; static char sol[64], best[64];
#define HB 26
static uint64_t *hk; static unsigned char *hg;
static int lb(const int *p){ int h=0;
  for(int i=0;i<K;i++) for(int j=i+1;j<K;j++){ int v=pair[i][j][p[i]][p[j]]; if(v>h) h=v; }
  int s=0; for(int d=0;d<10;d++){ int m=0; for(int i=0;i<K;i++) if(cnt[i][p[i]][d]>m) m=cnt[i][p[i]][d]; s+=m; }
  return s>h?s:h; }
static uint64_t key(const int *p){ uint64_t k=0; for(int i=0;i<K;i++) k=k*17+p[i]; return k; }
static int dfs(int *p, int g){
  nodes++;
  int done=1; for(int i=0;i<K;i++) if(p[i]<len[i]){done=0;break;}
  if(done){ memcpy(best,sol,g); best[g]=0; return 1; }
  if(g+lb(p)>L) return 0;
  uint64_t k=key(p); uint64_t h=(k*0x9E3779B97F4A7C15ULL)>>(64-HB);
  for(int t=0;t<8;t++){ uint64_t s=(h+t)&((1ULL<<HB)-1); if(hk[s]==k+1){ if(hg[s]<=g) return 0; hg[s]=g; goto go; } if(hk[s]==0){ hk[s]=k+1; hg[s]=g; goto go; } }
  go:;
  // candidate digits: next chars; order by how many words advance
  int c[10]={0}; for(int i=0;i<K;i++) if(p[i]<len[i]) c[w[i][p[i]]-'0']++;
  int ord[10], n=0; for(int d=0;d<10;d++) if(c[d]) ord[n++]=d;
  for(int a=0;a<n;a++) for(int b=a+1;b<n;b++) if(c[ord[b]]>c[ord[a]]){int t=ord[a];ord[a]=ord[b];ord[b]=t;}
  for(int t=0;t<n;t++){ int d=ord[t]; int q[K]; for(int i=0;i<K;i++) q[i]= (p[i]<len[i] && w[i][p[i]]-'0'==d)? p[i]+1 : p[i];
    sol[g]='0'+d; if(dfs(q,g+1)) return 1; }
  return 0; }
int main(int argc,char**argv){ L=atoi(argv[1]);
  const char *W[K]={"4242064951","5519385407","6886477529","3816749731","4700006089","9650560979","4371455527","5419320281","8761823531","3454277693"};
  for(int i=0;i<K;i++){ strcpy(w[i],W[i]); len[i]=strlen(W[i]); }
  for(int i=0;i<K;i++) for(int a=len[i];a>=0;a--) for(int d=0;d<10;d++) cnt[i][a][d]= a==len[i]?0: cnt[i][a+1][d] + (w[i][a]-'0'==d);
  for(int i=0;i<K;i++) for(int j=0;j<K;j++){ for(int a=len[i];a>=0;a--) for(int b=len[j];b>=0;b--){ int v;
      if(a==len[i]) v=len[j]-b; else if(b==len[j]) v=len[i]-a; else if(w[i][a]==w[j][b]) v=1+pair[i][j][a+1][b+1]; else { int x=pair[i][j][a+1][b], y=pair[i][j][a][b+1]; v=1+(x<y?x:y); }
      pair[i][j][a][b]=v; } }
  hk=calloc(1ULL<<HB,8); hg=malloc(1ULL<<HB);
  int p[K]={0}; printf("root lower bound %d\n", lb(p));
  int r=dfs(p,0);
  printf("L=%d: %s, nodes %lld\n", L, r?"FOUND":"NONE (exhaustive)", nodes); if(r) printf("%s\n",best);
  return 0; }
