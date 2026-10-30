// Puzzle 1139, general quadratic: v(n) = |a n^2 + b n + c|, halved when n is even (kept as 2 if it is 2).
// Usage: ./1139_quadratic a bmax cmin cmax R L0 SIEVE_LIMIT  -- prints every run of length >= L0 in n = -R..R.
// gcc -O2 -o 1139_quadratic 1139_quadratic.c
// Puzzle 1139 search: v(n) = |a n^2 + b n + c|, halved when n even.
// Parity analysis: valid iff a odd, b even, c = 2 mod 4 (then halved values are odd).
// For each (a,b,c) find the longest run of n in [-R,R] with all v(n) prime; report long runs.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static uint8_t *bits; static uint64_t LIM;
static int isp_small(uint64_t v){ if(v<2) return 0; if(v==2) return 1; if(!(v&1)) return 0; return (bits[v>>4]>>((v>>1)&7))&1; }
static uint64_t mulmod(uint64_t a,uint64_t b,uint64_t m){return (unsigned __int128)a*b%m;}
static uint64_t pw(uint64_t a,uint64_t e,uint64_t m){uint64_t r=1;a%=m;while(e){if(e&1)r=mulmod(r,a,m);a=mulmod(a,a,m);e>>=1;}return r;}
static int mr(uint64_t n){ if(n<LIM) return isp_small(n); if(!(n&1)) return 0;
  uint64_t d=n-1;int s=0;while(!(d&1)){d>>=1;s++;}
  static const uint64_t B[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){uint64_t a=B[i]%n; if(!a) continue; uint64_t x=pw(a,d,n); if(x==1||x==n-1) continue; int ok=0;
    for(int r=1;r<s;r++){x=mulmod(x,x,n); if(x==n-1){ok=1;break;}} if(!ok) return 0;} return 1; }
static long a_,b_,c_;
static uint64_t val(long n){ __int128 v=(__int128)a_*n*n+(__int128)b_*n+c_; if(v<0)v=-v; uint64_t u=(uint64_t)v; if(!(n&1) && u!=2) u>>=1; return u; }
static int ok(long n){ return mr(val(n)); }
int main(int argc,char**argv){
  long a=atol(argv[1]), bmax=atol(argv[2]); long cmin=atol(argv[3]), cmax=atol(argv[4]); long R=atol(argv[5]); int L0=atoi(argv[6]);
  LIM=(uint64_t)atoll(argv[7]);
  bits=calloc(LIM/16+1,1); // odd-only sieve
  memset(bits,0xff,LIM/16+1);
  for(uint64_t i=3;i*i<LIM;i+=2) if((bits[i>>4]>>((i>>1)&7))&1) for(uint64_t j=i*i;j<LIM;j+=2*i) bits[j>>4]&=~(1<<((j>>1)&7));
  bits[0]&=~1; // 1 not prime
  a_=a;
  for(long b=-bmax;b<=bmax;b+=2){ b_=b;
   long c0=cmin - ((cmin%4+4)%4) + 2;
   for(long c=c0;c<=cmax;c+=4){ c_=c;
    long lastend=-R-1;
    for(long an=-R+((L0-((-R)%L0+L0)%L0)%L0); an<=R; an+=L0){
      if(an<=lastend) continue;
      if(!ok(an)) continue;
      long lo=an,hi=an; while(lo-1>=-R && ok(lo-1)) lo--; while(hi+1<=R && ok(hi+1)) hi++;
      lastend=hi; int len=hi-lo+1;
      if(len>=L0){
        // distinct count
        static uint64_t vs[100000]; int m=0; for(long n=lo;n<=hi;n++) vs[m++]=val(n);
        int d=0; for(int i=0;i<m;i++){int dup=0; for(int j=0;j<i;j++) if(vs[j]==vs[i]){dup=1;break;} d+=!dup;}
        printf("a=%ld b=%ld c=%ld n=%ld..%ld len=%d distinct=%d\n",a,b,c,lo,hi,len,d); fflush(stdout);
      }
    }
   }
  }
  return 0;
}
