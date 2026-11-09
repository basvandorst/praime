// Puzzle 876. Build: gcc -O2 -o p876_runs p876_runs.c ; run: ./p876_runs N Q  (longest run of nonzero generated terms mod Q)
// max run of nonzero terms of t_m = t_{m-1}+...+t_{m-n} (mod q) over start vectors with at most one 0 (cap 400)
#include <stdio.h>
#include <stdlib.h>
int main(int argc,char**argv){ int n=atoi(argv[1]), q=atoi(argv[2]); int CAP=400;
  long long total=1; for(int i=0;i<n;i++) total*=q;
  int best=0; long long bestv=-1; int w[1000];
  for(long long v=0; v<total; v++){
    long long x=v; int z=0; for(int i=0;i<n;i++){ w[i]=x%q; x/=q; if(!w[i]) z++; } if(z>1) continue;
    int s=0; for(int i=0;i<n;i++) s+=w[i]; s%=q;
    int r=0, len=n;
    while(r<CAP){ int t=s; if(t==0) break; w[len]=t; s=(s+t-w[len-n]+q)%q; len++; r++; }
    if(r>best){best=r; bestv=v;}
    if(best>=CAP) break;
  }
  printf("n=%d q=%d: max run %d%s\n",n,q,best,best>=CAP?" (unbounded: a zero-free cycle exists)":"");
  return 0; }
