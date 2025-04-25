// minimal_standard.cpp
#include <iostream>
#include <ctime>
#include <iomanip>
using namespace std;

// Park–Miller constants and shuffle parameters
#define IA    16807
#define IM  2147483647
#define AM (1.0f/IM)
#define IQ    127773       // IM/IA
#define IR     2836       // IM%IA
#define NTAB      32
#define NDIV  (1 + (IM-1)/NTAB)
#define EPS   1.2e-7f
#define RNMX  (1.0f - EPS)

// ran1: Park–Miller minimal standard with Bays–Durham shuffle
float ran1(long *idum) {
    static long iy = 0, iv[NTAB];
    int j; long k;
    float temp;

    if (*idum <= 0 || iy == 0) {      // Initialize
        if (-(*idum) < 1) *idum = 1;  // avoid zero seed
        else *idum = -(*idum);
        for (j = NTAB+7; j >= 0; j--) {
            k = (*idum)/IQ;
            *idum = IA*(*idum - k*IQ) - IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ;                   // Update seed
    *idum = IA*(*idum - k*IQ) - IR*k;
    if (*idum < 0) *idum += IM;
    j = iy/NDIV;                      // Pick from shuffle table
    iy = iv[j];
    iv[j] = *idum;
    temp = AM * iy;                   // Scale to (0,1)
    return (temp > RNMX ? RNMX : temp);
}

int main() {
    long seed = -static_cast<long>(time(nullptr));
    cout << fixed << setprecision(6)
         << "Five uniform deviates from ran1:\n";
    for (int i = 0; i < 5; ++i)
        cout << "  " << ran1(&seed) << "\n";
    return 0;
}