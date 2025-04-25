#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;

//----------------------------------------------------------------------
// 1) ran0: Park–Miller “minimal standard” LCG (no shuffle)
//----------------------------------------------------------------------
#define IA0 16807
#define IM0 2147483647
#define IQ0 127773
#define IR0 2836
#define AM0 (1.0f/IM0)
#define MASK0 123459876

float ran0(long *idum) {
    long k = (*idum) / IQ0;
    *idum = IA0 * (*idum - k*IQ0) - IR0 * k;
    if (*idum < 0) *idum += IM0;
    return AM0 * (*idum);
}

//----------------------------------------------------------------------
// 2) ran1: Park–Miller with Bays–Durham shuffle
//----------------------------------------------------------------------
#define NTAB1 32
#define NDIV1 (1 + (IM0-1)/NTAB1)
#define EPS1 1.2e-7f
#define RNMX1 (1.0f - EPS1)

float ran1(long *idum) {
    static long iv[NTAB1], iy=0;
    if (*idum <= 0 || iy==0) {
        if (-(*idum) < 1) *idum = 1;
        else           *idum = -*idum;
        for (int j=NTAB1+7; j>=0; --j) {
            long k = (*idum)/IQ0;
            *idum = IA0*(*idum - k*IQ0) - IR0*k;
            if (*idum < 0) *idum += IM0;
            if (j<NTAB1) iv[j] = *idum;
        }
        iy = iv[0];
    }
    long k = (*idum)/IQ0;
    *idum = IA0*(*idum - k*IQ0) - IR0*k;
    if (*idum < 0) *idum += IM0;
    int j = iy/NDIV1;
    iy = iv[j];
    iv[j] = *idum;
    float temp = AM0 * iy;
    return (temp>RNMX1 ? RNMX1 : temp);
}

//----------------------------------------------------------------------
// 3) ran2: L’Ecuyer’s long‐period with Bays–Durham shuffle
//----------------------------------------------------------------------
#define IM1 2147483563
#define IM2 2147483399
#define AM2 (1.0f/IM1)
#define IMM12 (IM1-1)
#define IA1_2 40014
#define IA2_2 40692
#define IQ1_2 53668
#define IQ2_2 52774
#define IR1_2 12211
#define IR2_2 3791
#define NTAB2 32
#define NDIV2 (1 + IMM12/NTAB2)
#define EPS2 1.2e-7f
#define RNMX2 (1.0f - EPS2)

float ran2(long *idum) {
    static long idum2 = 123456789, iv[NTAB2], iy=0;
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum = 1;
        else              *idum = -*idum;
        idum2 = *idum;
        for (int j=NTAB2+7; j>=0; --j) {
            long k = (*idum)/IQ1_2;
            *idum = IA1_2*(*idum - k*IQ1_2) - IR1_2*k;
            if (*idum < 0) *idum += IM1;
            if (j<NTAB2) iv[j] = *idum;
        }
        iy = iv[0];
    }
    long k = (*idum)/IQ1_2;
    *idum = IA1_2*(*idum - k*IQ1_2) - IR1_2*k;
    if (*idum < 0) *idum += IM1;

    k = idum2/IQ2_2;
    idum2 = IA2_2*(idum2 - k*IQ2_2) - IR2_2*k;
    if (idum2 < 0) idum2 += IM2;

    int j = iy/NDIV2;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM12;
    float temp = AM2 * iy;
    return (temp>RNMX2 ? RNMX2 : temp);
}

//----------------------------------------------------------------------
// 4) ran3: Knuth’s subtractive generator
//----------------------------------------------------------------------
const int MBIG3 = 1000000000;
const int MZ3   = 0;
const float FAC3 = 1.0f/MBIG3;
static int inext3=0, inextp3=0;
static long ma3[56];
static bool init3 = false;

float ran3(long *idum) {
    if (*idum < 0 || !init3) {
        init3 = true;
        long mj = labs(161803398 - labs(*idum)) % MBIG3;
        ma3[55] = mj;
        long mk = 1;
        for (int i=1; i<=54; i++) {
            int ii = (21*i) % 55;
            ma3[ii] = mk;
            mk = mj - mk;
            if (mk < MZ3) mk += MBIG3;
            mj = ma3[ii];
        }
        for (int k=0; k<4; k++)
            for (int i=1; i<=55; i++) {
                ma3[i] -= ma3[1+(i+30)%55];
                if (ma3[i]<MZ3) ma3[i]+=MBIG3;
            }
        inext3 = 0;
        inextp3 = 31;
        *idum = 1;
    }
    if (++inext3 == 56) inext3 = 1;
    if (++inextp3== 56) inextp3= 1;
    long mj = ma3[inext3] - ma3[inextp3];
    if (mj < MZ3) mj += MBIG3;
    ma3[inext3] = mj;
    return mj * FAC3;
}

//----------------------------------------------------------------------
// Benchmark harness
//----------------------------------------------------------------------
int main(){
    const int N = 10000000;            // number of draws per generator
    volatile float sink;                 // prevent optimization

    // Seeds
    long seed0 = -static_cast<long>(time(nullptr));
    long seed1 = seed0 - 1;
    long seed2 = seed0 - 2;
    long seed3 = seed0 - 3;

    using Clock = chrono::high_resolution_clock;

    // Time ran0
    auto t0 = Clock::now();
    for(int i=0; i<N; ++i) sink = ran0(&seed0);
    float d0 = chrono::duration<float>(Clock::now() - t0).count();

    // Time ran1
    auto t1 = Clock::now();
    for(int i=0; i<N; ++i) sink = ran1(&seed1);
    float d1 = chrono::duration<float>(Clock::now() - t1).count();

    // Time ran2
    auto t2 = Clock::now();
    for(int i=0; i<N; ++i) sink = ran2(&seed2);
    float d2 = chrono::duration<float>(Clock::now() - t2).count();

    // Time ran3
    auto t3 = Clock::now();
    for(int i=0; i<N; ++i) sink = ran3(&seed3);
    float d3 = chrono::duration<float>(Clock::now() - t3).count();

    // Output relative timings
    cout << fixed << setprecision(3)
         << "ran0 = " << 1.0f          << "\n"
         << "ran1 = " << d1/d0         << "\n"
         << "ran2 = " << d2/d0         << "\n"
         << "ran3 = " << d3/d0         << "\n";

    return 0;
}
