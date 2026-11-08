// Puzzle 1279: primes p < B with ord_p(2) = d for the divisors d > 64 of L.
// Build: gcc -O2 -o p1279_primes p1279_primes.c ; run: ./p1279_primes 10080 4294967296 > big.txt
// Then: python3 p1279_primes_list.py 10080 big.txt > ../data/p1279_primes.txt (adds the orders d <= 64)
// list primes p < 2^32 with ord_p(2) = d for every divisor d > 64 of L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t pw(uint64_t b, uint64_t e, uint64_t m){uint64_t r=1;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
static int isprime(uint64_t n){ if(n<2)return 0; static const uint64_t a[]={2,3,5,7,11,13,17,19,23,29,31,37}; for(int i=0;i<12;i++){if(n%a[i]==0)return n==a[i];}
 uint64_t d=n-1;int s=0;while(!(d&1)){d>>=1;s++;} for(int i=0;i<12;i++){uint64_t x=pw(a[i],d,n); if(x==1||x==n-1)continue; int ok=0; for(int r=1;r<s;r++){x=(u128)x*x%n; if(x==n-1){ok=1;break;}} if(!ok)return 0;} return 1;}
int main(int argc,char**argv){ uint64_t L=strtoull(argv[1],0,10), B=strtoull(argv[2],0,10);
 int qs[16],nq=0; uint64_t t=L; for(uint64_t q=2;q*q<=t;q++) if(t%q==0){qs[nq++]=q; while(t%q==0)t/=q;} if(t>1)qs[nq++]=t;
 for(uint64_t d=65; d<=L; d++){ if(L%d) continue;
   uint64_t step = (d%2)? 2*d : d;
   for(uint64_t p=step+1; p<B; p+=step){ if(pw(2,d,p)!=1) continue; int ex=1; for(int i=0;i<nq;i++) if(d%qs[i]==0 && pw(2,d/qs[i],p)==1){ex=0;break;} if(!ex) continue; if(!isprime(p)) continue; printf("%llu %llu\n",(unsigned long long)p,(unsigned long long)d);} }
 return 0;}
