// count d-digit M ending in 9 with M + rev(M) a palindrome (exhaustive, d <= 10)
#include <stdio.h>
#include <stdint.h>
typedef unsigned long long u64;
static u64 rev(u64 x){u64 r=0;while(x){r=r*10+x%10;x/=10;}return r;}
static int pal(u64 x){return x==rev(x);}
int main(){
  u64 lo=1;
  for(int d=2; d<=10; d++){
    lo*=10; u64 L=lo/10, H=lo; long c=0; u64 ex=0;
    for(u64 M=L+9; M<H; M+=10){ if(pal(M+rev(M))){ c++; if(!ex) ex=M; } }
    printf("d=%2d: %ld M ending in 9 with M+rev(M) palindrome, first %llu\n", d, c, ex);
  }
  return 0;
}
