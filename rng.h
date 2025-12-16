// rng.h
#ifndef SI2530PROJECT_RNG_H
#define SI2530PROJECT_RNG_H

#include <random>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <thread>
#include <array>

namespace RNG {

    inline std::seed_seq make_seed_seq() {
        std::random_device rd;
        const auto tid = static_cast<unsigned long long>(
            std::hash<std::thread::id>{}(std::this_thread::get_id()));
        std::array<unsigned int, 4> seeds{
            rd(), rd(),
            static_cast<unsigned int>(tid & 0xFFFFFFFFu),
            static_cast<unsigned int>((tid >> 32) & 0xFFFFFFFFu)
        };
        return std::seed_seq(seeds.begin(), seeds.end());
    }

    // One independent engine per thread
    inline std::mt19937& engine() {
        // Make seed_seq an lvalue so mt19937 can bind to it
        static thread_local std::seed_seq seq = make_seed_seq();
        static thread_local std::mt19937 gen(seq);  // OK: lvalue seed_seq
        return gen;
    }

    inline std::size_t uniformIndex(std::size_t lo, std::size_t hi) {
        if (lo > hi) throw std::invalid_argument("uniformIndex: lo>hi");
        std::uniform_int_distribution<std::size_t> dist(lo, hi);
        return dist(engine());
    }

    inline int uniformSpin(int q) {
        if (q < 2) throw std::invalid_argument("uniformSpin: q<2");
        std::uniform_int_distribution<int> dist(0, q - 1);
        return dist(engine());
    }

    inline double uniformProb() {
        static thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(engine());
    }

    inline int uniformSpin(const std::size_t q) {
        if (q < 2 || q > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::invalid_argument("q out of range for uniformSpin");
        return uniformSpin(static_cast<int>(q));
    }
}

#endif // SI2530PROJECT_RNG_H
