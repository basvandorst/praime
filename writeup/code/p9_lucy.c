/* Puzzle 9: pi(x) and the sum of the primes <= x with Lucy's recursion, O(x^(3/4)) time, O(sqrt x) memory.
   Sums are exact in unsigned __int128.
   gcc -O2 -o p9_lucy writeup/code/p9_lucy.c -lm && ./p9_lucy 3531577135439   (~4 s) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static void p128(u128 x){char b[50];int i=49;b[i]=0;if(!x){puts("0");return;}while(x){b[--i]='0'+(int)(x%10);x/=10;}printf("%s",b+i);}
int main(int c,char**v){
  u64 n=strtoull(v[1],0,10); u64 r=(u64)sqrtl((long double)n); while(r*r>n)r--; while((r+1)*(r+1)<=n)r++;
  u64 m=r+n/r+5; u64*V=malloc(m*8); u64 cnt=0;
  for(u64 i=1;i<=r;i++) V[cnt++]=n/i;
  for(u64 x=V[cnt-1]-1;x>0;x--) V[cnt++]=x;
  u64 *C=malloc(cnt*8); u128 *S=malloc(cnt*16);
  for(u64 i=0;i<cnt;i++){u64 x=V[i]; C[i]=x-1; S[i]=(u128)x*(x+1)/2-1;}
  #define IDX(x) ((x)<=r ? cnt-(x) : n/(x)-1)
  for(u64 p=2;p<=r;p++){
    if(C[IDX(p)]==C[IDX(p-1)]) continue;
    u64 cp=C[IDX(p-1)]; u128 sp=S[IDX(p-1)]; u64 p2=p*p;
    for(u64 i=0;i<cnt && V[i]>=p2;i++){u64 j=IDX(V[i]/p); C[i]-=C[j]-cp; S[i]-=(u128)p*(S[j]-sp);}
  }
  printf("pi=%llu sum=",C[0]); p128(S[0]); puts("");
  return 0;
}
