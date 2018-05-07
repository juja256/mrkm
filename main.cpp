#include "primes.h"
#include "long_ar.h"
#include <iostream>
#include <time.h>

int main() {
    PrimeNumbers ctx(int(2048*2048*0.1));
    L_NUMBER p;
    GenStd prng(time(NULL));
    l_init_by_len(&p, 1024);

    ctx.FastPrimeMaurer( 1024, &p, &prng );
    std::cout << "Generated Prime: ";
    l_dump(&p, 'h');

    std::cout << "\nMR Check: " << ctx.MillerRabineTest(p, 10);
    
    return 0;
}
