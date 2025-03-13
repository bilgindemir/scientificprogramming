#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cstdlib>
#include <ctime>

namespace rangen {
    // Initialize the random seed. Call this once in main.
    inline void init_random() {
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
    }

    // Return a random float between a and b
    inline float frand(float a, float b) {
        return a + static_cast<float>(std::rand()) / (static_cast<float>(RAND_MAX/(b - a)));
    }
}

#endif // RANDOM_HPP
