// main.cpp
// assignment1_week6
// Created by Bilgin Demir on 4/15/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
using namespace std;

// chebft: compute n Chebyshev coefficients in c[0..n-1] for func over [a,b]
void chebft(float a, float b, vector<float>& c, int n, function<float(float)> func) {
    vector<float> f(n);
    float bma = 0.5f*(b - a), bpa = 0.5f*(b + a);
    for (int k = 0; k < n; k++) {
        float y = cosf(M_PI*(k + 0.5f)/n);
        f[k] = func(y*bma + bpa);
    }
    float fac = 2.0f/n;
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < n; k++)
            sum += f[k]*cos(M_PI*j*(k + 0.5f)/n);
        c[j] = fac * float(sum);
    }
}

// chebev: evaluate m-term Chebyshev series in c[0..m-1] at x over [a,b]
float chebev(float a, float b, const vector<float>& c, int m, float x) {
    if ((x - a)*(x - b) > 0.0f) 
        cout << "x not in range in chebev\n";
    float y  = (2.0f*x - a - b)/(b - a);
    float y2 = 2.0f*y;
    float d = 0.0f, dd = 0.0f, sv;
    for (int j = m-1; j >= 1; j--) {
        sv = d;
        d  = y2*d - dd + c[j];
        dd = sv;
    }
    return y*d - dd + 0.5f*c[0];
}

// Test function: e.g. exp(x)
float test_func(float x) {
    return expf(x);
}

int main() {
    float a = -1.0f, b = 1.0f;
    int n = 20;  // number of Chebyshev coefficients
    vector<float> c(n);

    // compute coefficients via chebft
    chebft(a, b, c, n, [](float x){ return test_func(x); });

    // print them
    cout << fixed << setprecision(6);
    cout << "Chebyshev coefficients:\n";
    for (int j = 0; j < n; j++) {
        cout << "c[" << j << "] = " << c[j] << "\n";
    }

    // test at a few points
    cout << "\nTesting approximation:\n";
    for (float x = a; x <= b; x += 0.5f) {
        float f_true = test_func(x);
        float f_cheb = chebev(a, b, c, n, x);
        cout << "x = " << setw(6) << x
             << "   true = " << setw(8) << f_true
             << "   cheb = " << setw(8) << f_cheb
             << "   err = " << setw(10) << f_cheb - f_true
             << "\n";
    }
    return 0;
}
