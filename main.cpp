#include "primes.h"
#include "long_ar.h"
#include <iostream>
#include <time.h>

int main() {
    PrimeNumbers ctx(int(2048*2048*0.1));
    L_NUMBER p;
    GenStd prng(time(NULL));
    l_init_by_len(&p, 2048);

    ctx.FastPrimeMaurer( 2048, &p, &prng );
    
    return 0;
}
