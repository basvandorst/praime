// Conjecture 111: local density of the shielded-primes sieve. Build: gcc -O2 -o conj111 conj-111.c -lm
// Run: ./conj111 X K 0   (K consecutive primes from X)   or   ./conj111 X K 4   (mod-6 random model from X)
// Local density of the shielded-primes sieve on K consecutive primes starting at X (segmented sieve),
// or (mode 4) on the mod-6 random model started at X. Ties left-to-right.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
static double *W; static int cmp_idx(const void*a,const void*b){ long i=*(const long*)a, j=*(const long*)b; double x=W[i], y=W[j]; return (x<y)?-1:(x>y)?1:((i<j)?-1:1); }
int main(int argc,char**argv){
  int64_t X=atoll(argv[1]); long K=atol(argv[2]); int mode=atoi(argv[3]);
  int64_t *P=malloc(sizeof(int64_t)*K); long c=0;
  if(mode==0){
    int64_t R=(int64_t)sqrt((double)X*1.01)+10; char *bs=calloc(R+1,1); int64_t *bp=malloc(sizeof(int64_t)*R); long nb=0;
    for(int64_t i=2;i<=R;i++) if(!bs[i]){ bp[nb++]=i; for(int64_t j=i*i;j<=R;j+=i) bs[j]=1; }
    int64_t S=1<<24; char *seg=malloc(S);
    for(int64_t lo=X; c<K; lo+=S){
      memset(seg,0,S);
      for(long t=0;t<nb;t++){ int64_t p=bp[t]; if(p*p>=lo+S) break; int64_t st=((lo+p-1)/p)*p; if(st<p*p) st=p*p; for(int64_t j=st;j<lo+S;j+=p) seg[j-lo]=1; }
      for(int64_t j=0;j<S && c<K;j++) if(!seg[j] && lo+j>1) P[c++]=lo+j;
    }
  } else { srand48(7); for(int64_t n=X|1; c<K; n+=2){ if(n%3==0) continue; if(drand48()<3.0/log((double)n)) P[c++]=n; } }
  W=malloc(sizeof(double)*K); long *idx=malloc(sizeof(long)*K);
  for(long i=0;i<K-1;i++){ idx[i]=i; W[i]=(double)(P[i+1]-P[i]); }
  qsort(idx,K-1,sizeof(long),cmp_idx);
  unsigned char *dead=calloc(K,1);
  for(long t=0;t<K-1;t++){ long i=idx[t]; if(!dead[i]&&!dead[i+1]) dead[i]=dead[i+1]=1; }
  long s=0; for(long i=0;i<K;i++) if(!dead[i]) s++;
  printf("X=%lld K=%ld mode=%d local density=%.7f  excess over e^-2 = %+.7f\n",(long long)X,K,mode,(double)s/K,(double)s/K-exp(-2.0)); return 0;
}
