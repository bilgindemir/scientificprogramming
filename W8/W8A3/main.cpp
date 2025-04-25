// ran2_demo.cpp
#include <iostream>
#include <ctime>
#include <cmath>
#include <iomanip>
using namespace std;

// Parameters from Numerical Recipes
#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0f/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1 + IMM1/NTAB)
#define EPS   1.2e-7f
#define RNMX  (1.0f - EPS)

// ran2: L’Ecuyer’s generator with Bays–Durham shuffle
float ran2(long *idum) {
    static long idum2 = 123456789;
    static long iy = 0, iv[NTAB];
    long j, k;
    float temp;

    if (*idum <= 0) {
        // Initialize shuffle table
        if (-(*idum) < 1) *idum = 1;
        else *idum = -*idum;
        idum2 = *idum;
        for (j = NTAB + 7; j >= 0; --j) {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k*IQ1) - IR1*k;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    // Update idum
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k*IQ1) - IR1*k;
    if (*idum < 0) *idum += IM1;

    // Update idum2
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k*IQ2) - IR2*k;
    if (idum2 < 0) idum2 += IM2;

    // Shuffle and form output
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;

    temp = AM * iy;
    return (temp > RNMX ? RNMX : temp);
}

int main(){
    // Seed with negative clock time to initialize
    long seed = -static_cast<long>(time(nullptr));

    cout << fixed << setprecision(6)
         << "Five uniform deviates from ran2:\n";
    for(int i = 0; i < 5; ++i)
        cout << "  " << ran2(&seed) << "\n";

    return 0;
}
