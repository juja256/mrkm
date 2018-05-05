#include "primes.h"
#include <time.h>
#include <math.h>
#include <algorithm>

static WORD base_pascal_seq_3[] = {1};
static WORD base_pascal_seq_5[] = {1};

#if ARCH == 64
static WORD base_pascal_seq_7[] = {1, 4, 2};
static WORD base_pascal_seq_11[] = {1, 4, 5, 9, 3};
static WORD base_pascal_seq_13[] = {1, 9, 3};
#else
static WORD base_pascal_seq_7[] = {1, 2, 4};
static WORD base_pascal_seq_11[] = {1, 9, 4, 3, 5};
static WORD base_pascal_seq_13[] = {1, 3, 9};
#endif


PrimeNumbers::PrimeNumbers(int n) {
    CalculatePrimes(n);
}

PrimeNumbers::PrimeNumbers() {}

void PrimeNumbers::CalculatePrimes(int n) {
    if (!primesTable.empty() && (n <= primesTable.back())) {
        return;
    }
    else if (!primesTable.empty() && (n > primesTable.back())) {
        bool isPrime = true;
        for (int i=primesTable.back()+1; i<=n; i++) {
            for (int j=2; j<=int(sqrt(i)); j++) {
                if (i % j == 0) {
                    isPrime = false;
                    break;
                }
            }
            if (isPrime) primesTable.push_back(i);
        }
    }
    else {
        bool isPrime = true;
        for (int i=2; i<=n; i++) {
            for (int j=2; j<=int(sqrt(i)); j++) {
                if (i % j == 0) {
                    isPrime = false;
                    break;
                }
            }
            if (isPrime) primesTable.push_back(i);
        }
    }
}

int PrimeNumbers::GetSmallPrimeNumber(int idx) {
    return primesTable.at(idx);
}

bool PrimeNumbers::IsSmallPrime(int n) {
    if (n<2) return false;
    for (int j=2; j<=int(sqrt(n)); j++) {
        if (n % j == 0) {
            return false;
        }
    }
    return true;
}

 /* OESIS A000720 */
int PrimeNumbers::GetSmallPrimesCount(int n) {
    static int small_pre[29] = {0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10};
    if (n <= 29) return small_pre[n];
    for (int i=n;i>0;i--) {
        if (IsSmallPrime(i)) {
            return std::distance(primesTable.begin(), std::find(primesTable.begin(), primesTable.end(), i));
        }
    }
}

static bool IsDivisible(L_NUMBER n, int p) {
    L_NUMBER tmp = {0,0};
    L_NUMBER r = {0,0};
    L_NUMBER zero = {0,0};
    l_copy(&r, &n);
    l_init(&tmp, n.len); tmp.words[0] = p;

    int k = l_bit_len(&n);
    int t = int( log2(p) ) + 1;

    if (k < t) {
        return false;
    }

    k = k-t;
    l_shift_l(&tmp, k, &tmp);

    while (k >= 0) {
        if ( l_sub(&r, &tmp, &r) == 1 ) {
            l_add(&r, &tmp, &r);
        }
        
        l_shift_r(&tmp, 1, &tmp);
        k--;
    }

    l_null(&tmp);

    bool res = l_cmp(&r, &tmp) == 0;

    l_free(&r);
    l_free(&tmp);
    return res;
}

/* From https://pdfs.semanticscholar.org/3df0/e64897ebbed46a6d034196f9b3962ca45b07.pdf */
void PrimeNumbers::FastPrimeMaurer(int k, L_NUMBER* N, PRNG* prng) {
    if (k<20) {
        unsigned n;
        bool isPrime;
        do {
            isPrime = true;
            n = prng->generateInt(k);
            int p_cnt = GetSmallPrimesCount( int( sqrt(n) ) );
            for (int i=0; i <= p_cnt; i++) {
                if (n % primesTable[i] == 0) {
                    isPrime = false;
                    break;
                }       
            }
        } while(isPrime);
        N->words[0] = n;
        return;
    } 
    else {
        float c = 0.1;
        int m =  20;
        float B = c * k * k;
        float r;
        if (k > 2*m) {
            do {
                float s = prng->generateFloat();
                r = pow(2, s-1);
            } while(int(k - r*k) > m);
        }
        else {
            r = 0.5;
        }
        FastPrimeMaurer(int(r*k), N, prng);
        bool success = false;
        L_NUMBER D, R, I;
        L_NUMBER unity;
        l_init(&unity, N->len); unity.words[0] = 1;
        l_init(&D, N->len); D.words[0] = 1;
        l_init(&R, N->len); R.words[0] = 1;
        l_init(&I, N->len);
        if (k-2 >= 0) {
            l_shift_l(&D, k-2, &D);
            l_div(&D, N, &I, &R);
            l_null(&R);
            R.words[0] = 1;
        }

        while (!success) {
            bool mayBePrime;
            L_NUMBER b1 = {0,0};
            L_NUMBER b2 = {0,0};
            L_NUMBER dd = {0,0};
            do {
                l_copy(&b1, &I);
                l_copy(&b2, &I);
                l_add(&b1, &unity, &b1);
                l_shift_l(&b2, 1, &b2);
                prng->generateBoundedLargeInt(b1, b2, &R);
                l_shift_l(&R, 1, &b1);
                l_mul(&b1, N, &dd);
                l_copy(&b1, &dd);
                l_add(&b1, &unity, &b1); // n = 2*R*q + 1
                mayBePrime = true;
                int p_cnt = GetSmallPrimesCount( int(B) );
                for (int i=0; i<=p_cnt; i++) {  
                    if (IsDivisible(b1, primesTable[i])) {
                        mayBePrime = false;
                        break;
                    }
                }
                
            } while( mayBePrime );
            L_NUMBER n = {0,0};
            l_copy(&n, &b1);
            l_null(&b2);
            b2.words[0] = 2;
            l_sub(&b1, &b2, &b1);
            L_NUMBER a = {0,0};
            L_NUMBER b = {0,0};
            L_NUMBER mu = {0,0};
            L_NUMBER pp = {0,0};
            l_init(&pp, N->len);
            prng->generateBoundedLargeInt(b2, b1, &a);
            m_pre_barret(2*N->len, &n, &mu);
            l_sub(&n, &unity, &pp);
            m_pow(&a, &pp, &n, &mu, &a); // a = a^{n-1} mod n
            if (l_cmp(&a, &unity) == 0) {
                
            }

            l_free(&b1);
            l_free(&b2);
            l_free(&dd);
            l_free(&a);
            
        }
    }
}

static void GenerateRandomForMR(L_NUMBER bound, L_NUMBER* out) {
    WORD len = bound.len;
    for (u32 i=0; i<len; i++) {
        out->words[i] = rand();
    }
    out->words[len - 1] %= bound.words[len - 1];
    if (!(out->words[len - 1])) out->words[len - 1]++;
}

static bool IsMutualPrimes(L_NUMBER a, L_NUMBER b) {
    bool fl;
    L_NUMBER t;
    l_init(&t, a.len);
    m_gcd(&a, &b, &t);
    if (t.words[0] == 1) fl = true;
    else fl = false;
    l_free(&t);
    return fl;
}

static bool PreDivisionTest(L_NUMBER p) {
    if (!(p.words[0] & 1)) return false; 
    WORD sum3 = 0;
    WORD sum5 = 0;
    WORD sum11 = 0;
    WORD sum7 = 0;
    WORD sum13 = 0;
    HALF* q = (HALF*)(p.words);
    for (u32 i=0; i<p.len*2; i++) {
        sum3 += q[i];
        sum5 += q[i];
        sum7 += q[i]*base_pascal_seq_7[i%3];
        sum11 += q[i]*base_pascal_seq_11[i%5];
        sum13 += q[i]*base_pascal_seq_13[i%3];
    }
    if ((sum3 % 3 == 0) || (sum5 % 5 == 0) || (sum7 % 7 == 0) || (sum11 % 11 == 0) || (sum13 % 13 == 0)) return false;
    return true;
}

bool PrimeNumbers::MillerRabineTest(L_NUMBER p, int k) {
    if (!PreDivisionTest(p)) { return false; }
    bool fl = false;
    srand(time(NULL));
    u32 n=0;
    L_NUMBER unity;
    L_NUMBER unity_inv;
    L_NUMBER tmp;
    L_NUMBER m;
    L_NUMBER d={0,0};
    L_NUMBER a;
    l_init(&a, p.len);
    l_init(&m, 2*p.len);
    l_init(&tmp, p.len);
    l_init(&unity, p.len); unity.words[0] = 1;
    l_init(&unity_inv, p.len); l_sub(&p, &unity, &unity_inv);
    l_copy(&d, &p);
    l_sub(&d, &unity, &d);
    WORD s=0;
    while (!(d.words[0] & 1)) {
        l_shift_r(&d, 1, &d);
        s++;

    }
    m_pre_barret(p.len*2, &p, &m);
    for (u32 n=0; n<k; n++) {
        GenerateRandomForMR(p, &a);

        if (IsMutualPrimes(a, p)) {
            m_pow(&a, &d, &p, &m, &tmp);
            if ((l_cmp(&unity, &tmp) == 0) || (l_cmp(&unity_inv, &tmp) == 0)) { 
                fl = true;
                continue;
            } 
            for (u32 i=1; i<s; i++) {
                m_sqr(&tmp, &p, &m, &tmp);
                if (l_cmp(&unity, &tmp) == 0) {
                    fl = false;
                    break;
                }
                else if (l_cmp(&unity_inv, &tmp) == 0) {
                    fl = true;
                    break;
                }
            }
            if (fl) continue;
            else break;
        }
        else {
            fl = false;
            break;
        }
    }
    l_free(&unity); l_free(&unity_inv); l_free(&a); l_free(&tmp); l_free(&m); l_free(&d);
    return fl;
}