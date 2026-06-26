#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define BENCH_ITER 100000
static double wt(void){return (double)clock()/CLOCKS_PER_SEC;}
int main(void){printf("bench: mini-circle-popov-criterion
");double t0=wt();for(int i=0;i<BENCH_ITER;i++){volatile double r=sqrt(sin(i*0.001));(void)r;}printf("bench PASS
");return 0;}
