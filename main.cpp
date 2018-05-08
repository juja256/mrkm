#include "primes.h"
#include "long_ar.h"
#include <iostream>

#ifdef _WIN64
#include <intrin.h>
#pragma intrinsic(_umul128) 
#include <windows.h>
#else
#include <x86intrin.h>
#include <time.h>
double GetTickCount(void) 
{
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now))
    return 0;
  return now.tv_sec * 1000.0 + now.tv_nsec / 1000000.0;
}
#endif // _WIN64

int main() {
    int bit_len = 1024;
    PrimeNumbers ctx(int(2048*2048*0.1));
    L_NUMBER p, q;
    GenStd prng_std(time(NULL));
    GenBM prng_bm(time(NULL));
    l_init_by_len(&p, bit_len);
    l_init_by_len(&q, bit_len);
    double t1 = GetTickCount();
    ctx.FastPrimeMaurer( bit_len, &p, &prng_std );
    double t2 = GetTickCount();
    ctx.PrimeMillerRabine(bit_len, &q, &prng_std );
    double t3 = GetTickCount();
    if (ctx.MillerRabineTest(p, 10)) {
        std::cout << "Maurer FastPrime Generated Prime(" << t2-t1 << " ms): ";
        l_dump(&p, 'h');
    }
    else {
        std::cerr << "Maurer Algorithm Error!\n";
    }
    std::cout << "Miller-Rabine Naїve Search Generated Prime(" << t3-t2 << " ms): ";
    l_dump(&q, 'h');
    
    return 0;
}
