#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

//-----------------------------
// Tuning parameters
//-----------------------------
#define ITMAX 100      // maximum iterations for gser/gcf
#define EPS    3e-7f   // relative accuracy
#define FPMIN  1e-30f  // avoid zero in continued fraction

//-----------------------------
// ln Gamma(x) via Lanczos
//-----------------------------
float gammln(float xx) {
    static double cof[6] = {
        76.18009172947146,  -86.50532032941677,
        24.01409824083091,  -1.231739572450155,
         0.1208650973866179e-2,  -0.5395239384953e-5
    };
    double x = xx, y = xx;
    double tmp = x + 5.5;
    tmp -= (x + 0.5)*log(tmp);
    double ser = 1.000000000190015;
    for(int j=0; j<6; j++) ser += cof[j]/(++y);
    return float(-tmp + log(2.5066282746310005*ser/x));
}

//-----------------------------
// Series representation of P(a,x)
//-----------------------------
void gser(float *gamser, float a, float x, float *gln) {
    *gln = gammln(a);
    if(x <= 0.0f) {
        *gamser = 0.0f;
        return;
    }
    float sum = 1.0f/a;
    float del = sum;
    float ap = a;
    for(int n=1; n<=ITMAX; n++) {
        ap += 1.0f;
        del *= x/ap;
        sum += del;
        if(fabs(del) < fabs(sum)*EPS) {
            *gamser = sum * expf(-x + a*logf(x) - *gln);
            return;
        }
    }
    cerr << "gser: a too large, ITMAX too small\n";
    *gamser = sum * expf(-x + a*logf(x) - *gln);
}

//-----------------------------
// Continued fraction for Q(a,x)
//-----------------------------
void gcf(float *gammcf, float a, float x, float *gln) {
    *gln = gammln(a);
    float b = x + 1.0f - a;
    float c = 1.0f / FPMIN;
    float d = 1.0f / b;
    float h = d;
    for(int i=1; i<=ITMAX; i++){
        float an = -i*(i - a);
        b += 2.0f;
        d = an*d + b;
        if(fabs(d) < FPMIN) d = FPMIN;
        c = b + an/c;
        if(fabs(c) < FPMIN) c = FPMIN;
        d = 1.0f/d;
        float del = d*c;
        h *= del;
        if(fabs(del - 1.0f) < EPS) break;
    }
    *gammcf = h * expf(-x + a*logf(x) - *gln);
}

//-----------------------------
// Incomplete gamma functions
// P(a,x) and Q(a,x)=1-P(a,x)
//-----------------------------
float gammp(float a, float x) {
    if(x < 0.0f || a <= 0.0f) {
        cerr << "gammp: invalid arguments\n";
        return 0.0f;
    }
    if(x < a+1.0f) {
        float gamser, gln;
        gser(&gamser,a,x,&gln);
        return gamser;
    } else {
        float gammcf, gln;
        gcf(&gammcf,a,x,&gln);
        return 1.0f - gammcf;
    }
}

float gammq(float a, float x) {
    if(x < 0.0f || a <= 0.0f) {
        cerr << "gammq: invalid arguments\n";
        return 0.0f;
    }
    if(x < a+1.0f) {
        float gamser, gln;
        gser(&gamser,a,x,&gln);
        return 1.0f - gamser;
    } else {
        float gammcf, gln;
        gcf(&gammcf,a,x,&gln);
        return gammcf;
    }
}

//-----------------------------
// Error function and complement
//-----------------------------
float erfF(float x) {
    return x<0.0f
         ? -gammp(0.5f, x*x)
         :  gammp(0.5f, x*x);
}

float erfFc(float x) {
    return x<0.0f
         ?  1.0f + gammp(0.5f, x*x)
         :  gammq(0.5f, x*x);
}

//-----------------------------
// Demo
//-----------------------------
int main(){
    cout << fixed << setprecision(7)
         << "  x     erfF(x)   std::erf   erfcF(x)   std::erfc\n";
    for(float x=-2.0f; x<=2.0f; x+=1.0f){
        cout
          << setw(5)<<x<<"   "
          << setw(8)<<erfF(x)<<"   "
          << setw(8)<<std::erf(x)<<"   "
          << setw(8)<<erfFc(x)<<"   "
          << setw(8)<<std::erfc(x)<<"\n";
    }
    return 0;
}
