// Puzzles 1079/1085. Build: gcc -O2 -o p1079_stones p1079_stones.c
// Run: ./p1079_stones N FIRST   (FIRST = 0: prime 3 east of the 1, FIRST = 1: south-east; the maximum over both is the answer)
// Prime stepping stones (puzzle 1079/1085): a 1 and n stones of value 2; place 3, 5, 7, ... in order, each equal to the
// sum of its already placed king-neighbours. Stones are placed lazily: when prime q goes to cell c, new stones may be put
// on empty cells next to c that touch no earlier prime cell. Exhaustive DFS; prints the largest prime reachable.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define G 80
#define C0 (G/2)
static int val[G][G], isprime_cell[G][G];
static int primes[100], np; static int FIRST; static int best=0, NST; static long long nodes=0;
static int bestgrid[G][G];
static int nb(int r,int c){ int s=0; for(int dr=-1;dr<=1;dr++) for(int dc=-1;dc<=1;dc++) if(dr||dc) s+=val[r+dr][c+dc]; return s; }
static int touches_prime(int r,int c){ for(int dr=-1;dr<=1;dr++) for(int dc=-1;dc<=1;dc++) if((dr||dc)&&isprime_cell[r+dr][c+dc]) return 1; return 0; }
static int has_odd_nb(int r,int c){ for(int dr=-1;dr<=1;dr++) for(int dc=-1;dc<=1;dc++) if((dr||dc)&&(val[r+dr][c+dc]&1)) return 1; return 0; }
static int minr=C0,maxr=C0,minc=C0,maxc=C0;
static void dfs(int qi, int left){
  nodes++;
  if(primes[qi-1] > best && qi>0){ best=primes[qi-1]; memcpy(bestgrid,val,sizeof val); }
  int q=primes[qi];
  int r0=minr-1, r1=maxr+1, c0=minc-1, c1=maxc+1;
  if(r0<2||c0<2||r1>G-3||c1>G-3){ fprintf(stderr,"grid too small\n"); exit(1); }
  for(int r=r0;r<=r1;r++) for(int c=c0;c<=c1;c++){
    if(val[r][c]) continue;
    if(qi==0 && !(FIRST==0 ? (r==C0 && c==C0+1) : (r==C0+1 && c==C0+1))) continue;   // symmetry class of the first move
    if(!has_odd_nb(r,c)) continue;
    int s0=nb(r,c); if(s0>q) continue;
    int d=q-s0; if(d&1) continue; int k=d/2; if(k>left) continue;
    // candidate stone cells
    int cr[8],cc[8],m=0;
    for(int dr=-1;dr<=1;dr++) for(int dc=-1;dc<=1;dc++) if(dr||dc){ int a=r+dr,b=c+dc; if(!val[a][b] && !touches_prime(a,b)) { cr[m]=a; cc[m]=b; m++; } }
    if(k>m) continue;
    // enumerate subsets of size k
    for(int mask=0;mask<(1<<m);mask++){ if(__builtin_popcount(mask)!=k) continue;
      int sv[4]={minr,maxr,minc,maxc};
      for(int t=0;t<m;t++) if(mask>>t&1){ val[cr[t]][cc[t]]=2; }
      val[r][c]=q; isprime_cell[r][c]=1;
      for(int t=0;t<m;t++) if(mask>>t&1){ if(cr[t]<minr)minr=cr[t]; if(cr[t]>maxr)maxr=cr[t]; if(cc[t]<minc)minc=cc[t]; if(cc[t]>maxc)maxc=cc[t]; }
      if(r<minr)minr=r; if(r>maxr)maxr=r; if(c<minc)minc=c; if(c>maxc)maxc=c;
      dfs(qi+1,left-k);
      val[r][c]=0; isprime_cell[r][c]=0; for(int t=0;t<m;t++) if(mask>>t&1) val[cr[t]][cc[t]]=0;
      minr=sv[0];maxr=sv[1];minc=sv[2];maxc=sv[3];
    }
  }
}
int main(int argc,char**argv){ NST=atoi(argv[1]); FIRST=atoi(argv[2]);
  for(int p=3;p<500;p++){ int ok=1; for(int d=2;d*d<=p;d++) if(p%d==0){ok=0;break;} if(ok) primes[np++]=p; }
  memset(val,0,sizeof val); val[C0][C0]=1;
  dfs(0,NST);
  printf("n=%d: largest prime %d, nodes %lld\n",NST,best,nodes);
  int r0=G,r1=0,c0=G,c1=0; for(int r=0;r<G;r++) for(int c=0;c<G;c++) if(bestgrid[r][c]){ if(r<r0)r0=r; if(r>r1)r1=r; if(c<c0)c0=c; if(c>c1)c1=c; }
  for(int r=r0;r<=r1;r++){ for(int c=c0;c<=c1;c++){ if(bestgrid[r][c]) printf("%4d",bestgrid[r][c]); else printf("   ."); } printf("\n"); }
  return 0; }
