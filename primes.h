#ifndef PRIMES_H
#define PRIMES_H

#include "long_ar.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

typedef unsigned char BIT;
typedef unsigned char BYTE;

class PRNG {
public:
    virtual BIT next()=0;
    virtual void generateByteArray(BYTE* arr, int bit_size) {
        int l = (bit_size%8 == 0) ? (bit_size/8) : (bit_size/8 + 1);
        memset(arr, 0, l);
        for (int i=0; i<bit_size; i++) {
            arr[i/8] |= next() << (i%8);
        }
    }
    virtual unsigned generateInt(int bit_len) {
        if (bit_len > 32) {
            bit_len = 32;
        }
        unsigned res = 0;
        for (int i=0; i<bit_len; i++) {
            res |= next() << i;
        }
        return res;
    }
    virtual float generateFloat() {
        return generateInt(32) * 1.0 / (unsigned)(-1);
    }
    virtual void generateBoundedLargeInt(L_NUMBER b1, L_NUMBER b2, L_NUMBER* res) {
        L_NUMBER tmp;
        if (res->len == 0) l_init(res, b2.len);
        l_init(&tmp, b1.len);
        l_sub(&b2, &b1, &tmp);
        int l = l_bit_len(&tmp);
        generateByteArray((BYTE*)(res->words), l);
        if (l_cmp(res, &tmp) == 1) l_sub(res, &tmp, res);
        l_add(&b1, res, res);
        l_free(&tmp);
    }
};

class Gen55 : public PRNG {
    int s;

    public:
    Gen55() {
        s = 1;
    }

    BIT next() override {
        s ^= 1;
        return s;
    }
};

class GenStd : public PRNG {
    int s;
    public:
    GenStd(int s_) {
        srand(s_);
        s = rand();
    }
    BIT next() override {
        s = rand();
        return s & 1;
    }
};

class GenBM : public PRNG {
    L_NUMBER a, p, q, t, mu;
    public:

    GenBM() {
        l_init_by_str(&a, "0x5B88C41246790891C095E2878880342E88C79974303BD0400B090FE38A688356");
        l_init_by_str(&p, "0xCEA42B987C44FA642D80AD9F51F10457690DEF10C83D0BC1BCEE12FC3B6093E3");
        l_init_by_str(&q, "0x675215CC3E227D3216C056CFA8F8822BB486F788641E85E0DE77097E1DB049F1");
        l_init_by_len(&t, 256);
        t.words[0] = 1;
        l_init_by_len(&mu, 512);
	    m_pre_barret(8, &p, &mu);
    }

    GenBM(L_NUMBER seed) : GenBM() {
        setSeed(seed);
    }

    ~GenBM() {
        l_free(&a); l_free(&p); l_free(&q); l_free(&mu); l_free(&t);
    }

    void setSeed(L_NUMBER seed) {
        l_copy(&t, &seed);
        m_pow(&a, &t, &p, &mu, &t);
    }

    BIT next() override {
        m_pow(&a, &t, &p, &mu, &t);
        return (l_cmp(&q, &t) == 1) ? 1 : 0;
    }
};

class PrimeNumbers {
    std::vector<int> primesTable;
public:
    PrimeNumbers();
    PrimeNumbers(int n);
    bool IsSmallPrime(int n);
    void CalculatePrimes(int n);
    int GetSmallPrimeNumber(int idx);
    int GetSmallPrimesCount(int n);
    bool MillerRabineTest(L_NUMBER p, int k);
    void FastPrimeMaurer(int k, L_NUMBER* N, PRNG* prng);
};


#endif