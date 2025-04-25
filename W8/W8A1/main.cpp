// uniform_deviates.cpp
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;

// ----- 1) System-supplied rand() wrapper -----
inline float uni_sys() {
    // rand() returns [0 .. RAND_MAX], so divide by RAND_MAX+1
    return rand() / (RAND_MAX + 1.0f);
}

// ----- 2) Park–Miller minimal standard RNG (ran0) -----

// Constants for Schrage’s method: a=16807, m=2^31−1
#define IA 16807
#define IM 2147483647
#define IQ 127773          // IM / IA
#define IR 2836            // IM % IA
#define AM (1.0f/IM)
#define MASK 123459876     // to prevent idum=0 corner-case

float ran0(long *idum) {
    long k;
    float ans;
    *idum ^= MASK;                      // mix in MASK so seed=0 still works
    k = (*idum) / IQ;
    *idum = IA * (*idum - k*IQ) - IR * k; 
    if (*idum < 0) *idum += IM;         // ensure positive
    ans = AM * (*idum);                 // scale to (0,1)
    *idum ^= MASK;                      // undo mask
    return ans;
}

int main() {
    // seed both generators from the clock
    unsigned seed = static_cast<unsigned>(time(nullptr));
    srand(seed);

    long idum = -static_cast<long>(seed);  // ran0 wants negative seed to initialize

    cout << fixed << setprecision(6)
         << "Five uniform deviates from std::rand():\n";
    for (int i = 0; i < 5; ++i)
        cout << "  " << uni_sys() << "\n";

    cout << "\nFive uniform deviates from ran0 (Park–Miller):\n";
    for (int i = 0; i < 5; ++i)
        cout << "  " << ran0(&idum) << "\n";

    return 0;
}
