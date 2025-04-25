// ran3_knuth.cpp
#include <iostream>
#include <cmath>
#include <ctime>
#include <iomanip>
using namespace std;

// Constants per Knuth’s Algorithm
const int MBIG = 1000000000;
const int MSEED = 161803398;
const int MZ    = 0;
const float FAC = 1.0f/MBIG;

// ran3: Knuth’s subtractive generator
float ran3(long *idum) {
    static int inext = 0, inextp = 0;
    static long ma[56];
    static bool initialized = false;

    if (*idum < 0 || !initialized) {
        // Initialization
        initialized = true;
        long mj = labs(MSEED - labs(*idum)) % MBIG;
        ma[55] = mj;
        long mk = 1;
        for (int i = 1; i <= 54; i++) {
            int ii = (21*i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) mk += MBIG;
            mj = ma[ii];
        }
        // Warm up
        for (int k = 0; k < 4; k++)
            for (int i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i+30)%55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext = 0;
        inextp = 31;
        *idum = 1;  // seed is now “used”
    }

    // Generate a new random value
    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;

    long mj = ma[inext] - ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext] = mj;

    return mj * FAC;
}

int main() {
    // Seed with negative time to trigger initialization
    long seed = -static_cast<long>(time(nullptr));

    cout << fixed << setprecision(6)
         << "Five uniform deviates from ran3 (Knuth):\n";
    for (int i = 0; i < 5; ++i) {
        cout << "  " << ran3(&seed) << "\n";
    }

    return 0;
}
