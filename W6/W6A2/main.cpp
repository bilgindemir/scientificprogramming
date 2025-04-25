#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
using namespace std;

// --- Chebyshev fit (from Assignment 1) ---
void chebft(float a, float b, vector<float>& c, int n, function<float(float)> func) {
    vector<float> f(n);
    float bma = 0.5f*(b - a), bpa = 0.5f*(b + a);
    for (int k = 0; k < n; ++k) {
        float y = cosf(M_PI*(k + 0.5f)/n);
        f[k] = func(y*bma + bpa);
    }
    float fac = 2.0f/n;
    for (int j = 0; j < n; ++j) {
        double sum = 0.0;
        for (int k = 0; k < n; ++k)
            sum += f[k]*cos(M_PI*j*(k + 0.5f)/n);
        c[j] = fac * float(sum);
    }
}

float chebev(float a, float b, const vector<float>& c, int m, float x) {
    if ((x - a)*(x - b) > 0.0f)
        cout << "Warning: x out of range in chebev\n";
    float y  = (2.0f*x - a - b)/(b - a);
    float y2 = 2.0f*y;
    float d = 0.0f, dd = 0.0f, sv;
    for (int j = m-1; j >= 1; --j) {
        sv = d;
        d  = y2*d - dd + c[j];
        dd = sv;
    }
    return y*d - dd + 0.5f*c[0];
}

// --- Assignment 2: Derivative coefficients chder ---
void chder(float a, float b, const vector<float>& c, vector<float>& cder) {
    int n = c.size();
    cder.assign(n, 0.0f);
    if (n >= 2) {
        cder[n-1]   = 0.0f;
        cder[n-2] = 2.0f*(n-1)*c[n-1];
        for (int j = n-3; j >= 0; --j)
            cder[j] = cder[j+2] + 2.0f*(j+1)*c[j+1];
    }
    float con = 2.0f/(b - a);
    for (int j = 0; j < n; ++j)
        cder[j] *= con;
}

// --- Assignment 2: Integral coefficients chint ---
void chint(float a, float b, const vector<float>& c, vector<float>& cint) {
    int n = c.size();
    cint.assign(n, 0.0f);
    if (n >= 2) {
        float sum = 0.0f, fac = 1.0f;
        float con = 0.25f*(b - a);
        for (int j = 1; j <= n-2; ++j) {
            cint[j] = con*(c[j-1] - c[j+1])/float(j);
            sum    += fac*cint[j];
            fac     = -fac;
        }
        cint[n-1] = con*c[n-2]/float(n-1);
        sum      += fac*cint[n-1];
        cint[0]   = 2.0f*sum;
    }
}

// --- Test function and true derivatives/integrals ---
float f_true(float x)    { return expf(x); }
float df_true(float x)   { return expf(x); }
float Fint_true(float a, float x) { return expf(x) - expf(a); }

int main(){
    // Interval and polynomial order
    float a = -1.0f, b = 1.0f;
    int   n = 20;
    vector<float> c(n), cder, cint;

    // Compute Chebyshev coeffs for f(x)=e^x
    chebft(a, b, c, n, f_true);

    // Derivative and integral coefficients
    chder(a, b, c, cder);
    chint(a, b, c, cint);

    cout << fixed << setprecision(6);

    // Test at sample points
    cout << " x      f'≈      f' true     error    ∫f≈      ∫f true    error\n";
    for (float x = a; x <= b; x += 0.5f) {
        float df_approx = chebev(a, b, cder, n, x);
        float df_exp    = df_true(x);
        float Fi_approx = chebev(a, b, cint, n, x);
        float Fi_true   = Fint_true(a, x);
        cout
          << setw(6)<<x << "  "
          << setw(8)<<df_approx << "  "
          << setw(8)<<df_exp    << "  "
          << setw(8)<<(df_approx-df_exp) << "  "
          << setw(8)<<Fi_approx << "  "
          << setw(8)<<Fi_true   << "  "
          << setw(8)<<(Fi_approx-Fi_true)
          << "\n";
    }
    return 0;
}
