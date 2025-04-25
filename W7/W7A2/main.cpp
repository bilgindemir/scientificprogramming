#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

// Tuning parameters
const int    MAXIT = 100;      // max iterations
const float  EPS   = 1e-7f;    // relative accuracy
const float  FPMIN = 1e-30f;   // to avoid underflow
const float  EULER= 0.57721566490153286060f; // Euler’s constant

// expint: computes E_n(x)
float expint(int n, float x) {
    if (n < 0 || x < 0.0f || (x == 0.0f && (n==0 || n==1))) {
        cerr << "bad arguments in expint\n";
        return 0.0f;
    }
    int nm1 = n - 1;
    // Special cases
    if (n == 0) return expf(-x)/x;
    if (x == 0.0f) return 1.0f/nm1;

    // Continued fraction for x>1
    if (x > 1.0f) {
        float b = x + n;
        float c = 1.0f / FPMIN;
        float d = 1.0f / b;
        float h = d;
        for (int i = 1; i <= MAXIT; i++) {
            float a = -i*(nm1 + i);
            b += 2.0f;
            d = a*d + b;
            if (fabsf(d) < FPMIN) d = FPMIN;
            c = b + a/c;
            if (fabsf(c) < FPMIN) c = FPMIN;
            d = 1.0f / d;
            float del = c*d;
            h *= del;
            if (fabsf(del - 1.0f) < EPS) {
                return h * expf(-x);
            }
        }
        cerr << "continued fraction failed in expint\n";
        return h * expf(-x);
    }

    // Power-series for x<=1
    float ans = (nm1 != 0 ? 1.0f/nm1 : -logf(x) - EULER);
    float fact = 1.0f;
    for (int i = 1; i <= MAXIT; i++) {
        fact *= -x / i;
        float del;
        if (i != nm1) {
            del = -fact / (i - nm1);
        } else {
            // digamma for integer n
            float psi = -EULER;
            for (int ii = 1; ii <= nm1; ii++)
                psi += 1.0f/ii;
            del = fact * (-logf(x) + psi);
        }
        ans += del;
        if (fabsf(del) < fabsf(ans)*EPS)
            return ans;
    }
    cerr << "series failed in expint\n";
    return ans;
}

// ei: computes Ei(x) for x>0
float ei(float x) {
    if (x <= 0.0f) {
        cerr << "Bad argument in ei\n";
        return NAN;
    }
    if (x < FPMIN) {
        // avoid division-by-zero
        return logf(x) + EULER;
    }
    if (x <= -logf(EPS)) {
        // Power series for small x
        float sum = 0.0f, fact = 1.0f;
        for (int k = 1; k <= MAXIT; k++) {
            fact *= x/k;
            float term = fact/k;
            sum += term;
            if (term < EPS*sum)
                return sum + logf(x) + EULER;
        }
        cerr << "Series failed in ei\n";
        return sum + logf(x) + EULER;
    }
    // Asymptotic series for large x
    float sum = 0.0f, term = 1.0f, prev;
    for (int k = 1; k <= MAXIT; k++) {
        prev = term;
        term *= k/x;
        if (term < EPS) break;
        if (term < prev)
            sum += term;
        else {
            sum -= prev;
            break;
        }
    }
    return expf(x)*(1.0f + sum)/x;
}

int main(){
    cout << fixed << setprecision(8);

    // Demonstrate E_n(x) for n=1,2,3 at x=0.5 and x=2.0
    for (int n : {1,2,3}) {
        for (float x : {0.5f, 2.0f}) {
            cout << "E_"<<n<<"("<<x<<") = "<<setw(12)
                 <<expint(n,x)<<"\n";
        }
    }

    // Demonstrate Ei(x) for x=0.5 and x=5.0
    for (float x : {0.5f, 5.0f}) {
        cout << "Ei("<<x<<") = "<<setw(12)<<ei(x)<<"\n";
    }

    return 0;
}
