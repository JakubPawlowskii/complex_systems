//
// Created by Jakub on 23.02.2021.
// Class handling all-purpose pseudo-random number generation
// Header only
//

#ifndef SIMULATION_PRNG_H
#define SIMULATION_PRNG_H
#include <random>

class prng {
private:
    /**
    * Utilities for PRNG using xoshiro256**
    * http://prng.di.unimi.it/xoshiro256starstar.c
    */
    std::random_device device;
    uint64_t s[4]{};

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    uint64_t next() {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;

        const uint64_t t = s[1] << 17;

        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];

        s[2] ^= t;

        s[3] = rotl(s[3], 45);

        return result;
    }

public:
    virtual ~prng() = default;

    prng(){
        for(int i = 0; i < 3; i++)
            s[i] = (static_cast<uint64_t>(device()) << 32) | device();
    }

    /*
 * Return random 64-bit integer
 */
    uint64_t random_uint_64() {
        return next();
    }

/*
 * Return random double from [0,1)
 */
    double random_double() {
        uint64_t x = next();
        const union {
            uint64_t i;
            double d;
        } u = {UINT64_C(0x3FF) << 52 | x >> 12};
        return u.d - 1.0;
    }


};


#endif //SIMULATION_PRNG_H
