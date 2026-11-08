// Puzzle 771. Build: gcc -O2 -o p771_anneal p771_anneal.c -lm ; run: ./p771_anneal N P ITERATIONS SEED
// (N = 6, P = 293, 4e7 iterations finds a grid in seconds; N = 5, P = 227 always ends one prime short)
// Simulated annealing: N x N digit grid; every prime <= P must be readable along a straight run (8 directions).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static int N, P, nr; static int run[4000][4], rl[4000];
static int primes[400], np;
static unsigned long long rs=88172645463325252ULL; static unsigned rnd(){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (unsigned)(rs>>11); }
static int cost(const int *g){ static char rd[10000]; memset(rd,0,10000);
  for(int r=0;r<nr;r++){ int v=0; if(rl[r]>1 && g[run[r][0]]==0) continue; for(int t=0;t<rl[r];t++) v=v*10+g[run[r][t]]; rd[v]=1; }
  int c=0; for(int i=0;i<np;i++) if(!rd[primes[i]]) c++; return c; }
int main(int argc,char**argv){ N=atoi(argv[1]); P=atoi(argv[2]); long long iters=atoll(argv[3]); rs^=atoll(argv[4])*0x9E3779B97F4A7C15ULL;
  for(int p=2;p<=P;p++){ int ok=1; for(int d=2;d*d<=p;d++) if(p%d==0){ok=0;break;} if(ok) primes[np++]=p; }
  int maxlen = P>=100?3:(P>=10?2:1); if(P>=1000) maxlen=4;
  int dr[8]={0,0,1,-1,1,1,-1,-1}, dc[8]={1,-1,0,0,1,-1,1,-1}; nr=0;
  for(int r=0;r<N;r++) for(int c=0;c<N;c++){ run[nr][0]=r*N+c; rl[nr++]=1;
    for(int L=2;L<=maxlen;L++) for(int d=0;d<8;d++){ int ok=1; for(int t=0;t<L;t++){ int rr=r+dr[d]*t, cc=c+dc[d]*t; if(rr<0||rr>=N||cc<0||cc>=N){ok=0;break;} run[nr][t]=rr*N+cc; } if(ok) rl[nr++]=L; } }
  int g[64]; for(int i=0;i<N*N;i++) g[i]=rnd()%10;
  int cur=cost(g), best=cur; int bg[64]; memcpy(bg,g,sizeof g); double T=2.0;
  for(long long it=0;it<iters && best>0;it++){
    T = 2.0*pow(0.02/2.0,(double)it/iters);
    int i=rnd()%(N*N), old=g[i], nd=rnd()%10; if(nd==old) continue; g[i]=nd;
    int c=cost(g);
    if(c<=cur || exp((cur-c)/T) > (rnd()%100000)/100000.0){ cur=c; if(c<best){best=c; memcpy(bg,g,sizeof g);} } else g[i]=old;
  }
  printf("N=%d P=%d best missing=%d grid=",N,P,best); for(int r=0;r<N;r++){ for(int c=0;c<N;c++) printf("%d",bg[r*N+c]); printf(r<N-1?"/":"\n"); }
  return 0; }
