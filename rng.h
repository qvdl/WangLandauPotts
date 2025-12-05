//
// Created by quinten on 12/5/25.
//

#ifndef SI2530PROJECT_RNG_H
#define SI2530PROJECT_RNG_H


#include <random>
#include <cstddef>

namespace RNG {
    // The central engine used everywhere
    inline std::mt19937& engine() {
        static std::mt19937 gen(std::random_device{}());
        return gen;
    }

    // Uniform integer in [min, max]
    inline std::size_t uniformIndex(const std::size_t min, const std::size_t max) {
        std::uniform_int_distribution<std::size_t> dist(min, max);
        return dist(engine());
    }

    // Uniform integer in [0, q-1]
    inline int uniformSpin(const int q) {
        std::uniform_int_distribution<int> dist(0, q - 1);
        return dist(engine());
    }

    // Uniform floating-point in [0,1)
    inline double uniformProb() {
        static std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(engine());
    }
}




#endif //SI2530PROJECT_RNG_H