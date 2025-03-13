// main.cpp
// assignment3_week3
// Created by Bilgin Demir on 3/3/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>
#include <iomanip>

#define EPS 1.0e-6       // Desired fractional accuracy
#define JMAX 20          // Maximum iterations for qtrap and qsimp
#define JMAXP 20         // Maximum iterations for qromb
#define K 5              // Use 5 points for extrapolation in qromb

// Global counter to track how many times trapzd is called.
static int trapzdCount = 0;

// Integrand: f(x) = x^4 * log(x + sqrt(x^2 + 1))
// over the interval [0, 2].
float func(float x) {
    return x*x*x*x * std::log(x + std::sqrt(x*x + 1.0f));
}

// polint: Polynomial interpolation used in Romberg extrapolation.
// This version expects arrays with 1-indexing (i.e. elements stored in indices 1..n).
void polint(const float xa[], const float ya[], int n, float x, float &y, float &dy) {
    int i, m, ns = 1;  // start at index 1
    float den, dif, dift, ho, hp, w;
    float *c = new float[n+1];
    float *d = new float[n+1];

    dif = std::fabs(x - xa[1]);
    for (i = 1; i <= n; i++) {
        dift = std::fabs(x - xa[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns];
    ns--;
    for (m = 1; m < n; m++) {
        for (i = 1; i <= n - m; i++) {
            ho = xa[i]   - x;
            hp = xa[i+m] - x;
            w  = c[i+1]  - d[i];
            den = ho - hp;
            if (den == 0.0f) {
                std::cerr << "Error in polint: den=0\n";
                delete[] c; delete[] d;
                return;
            }
            den = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }
        if (2*ns < (n - m))
            dy = c[ns+1];
        else {
            dy = d[ns];
            ns--;
        }
        y += dy;
    }
    delete[] c;
    delete[] d;
}

// trapzd: One stage of refinement of the trapezoidal rule.
// EXACT replication of Numerical Recipes code using 1-indexed logic.
float trapzd(float (*f)(float), float a, float b, int n) {
    trapzdCount++;  // Increment counter

    static float s;
    int it, j;
    float x, tnm, sum, del;

    if (n == 1) {
        s = 0.5f * (b - a) * (f(a) + f(b));
        return s;
    } else {
        for (it = 1, j = 1; j < n - 1; j++) {
            it <<= 1;  // it = 2^(n-2)
        }
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5f * del;
        sum = 0.0f;
        for (j = 1; j <= it; j++, x += del)
            sum += f(x);
        s = 0.5f * (s + (b - a) * sum / tnm);
        return s;
    }
}

// qtrap: Basic trapezoidal integration driver.
// Uses trapzd repeatedly until convergence. Expected to require 13 calls.
float qtrap(float (*f)(float), float a, float b) {
    float s, olds = 0.0f;
    for (int j = 1; j <= JMAX; j++) {
        s = trapzd(f, a, b, j);
        if (j > 5) {
            if (std::fabs(s - olds) < EPS * std::fabs(olds))
                return s;
        }
        olds = s;
    }
    std::cerr << "Too many steps in qtrap\n";
    return s;
}

// qsimp: Simpson's rule integration driver.
// Uses trapzd plus a Richardson extrapolation step.
// Expected to require 8 calls.
float qsimp(float (*f)(float), float a, float b) {
    float s, st, ost = 0.0f, os = 0.0f;
    for (int j = 1; j <= JMAX; j++) {
        st = trapzd(f, a, b, j);
        s = (4.0f * st - ost) / 3.0f;
        if (j > 5) {
            if (std::fabs(s - os) < EPS * std::fabs(os))
                return s;
        }
        os = s;
        ost = st;
    }
    std::cerr << "Too many steps in qsimp\n";
    return 0.0f;
}

// qromb: Romberg integration driver.
// We want it to converge on the very first extrapolation: after 5 calls to trapzd.
// To do this, we use 1-indexed arrays for h and s (of size JMAXP+1 and JMAXP+2),
// and set K = 5 so that when j == 5 we attempt extrapolation using s[1..5] and h[1..5].
float qromb(float (*f)(float), float a, float b) {
    float ss, dss;
    // Allocate arrays with 1-indexing: ignore index 0.
    float s[JMAXP+1], h[JMAXP+2];
    h[1] = 1.0f;
    for (int j = 1; j <= JMAXP; j++) {
        s[j] = trapzd(f, a, b, j);
        if (j >= K) {
            // Pass pointers so that polint sees elements 1..K:
            polint(h + j - K + 0, s + j - K + 0, K, 0.0f, ss, dss);
            if (std::fabs(dss) <= EPS * std::fabs(ss))
                return ss;
        }
        h[j+1] = 0.25f * h[j];
    }
    std::cerr << "Too many steps in qromb\n";
    return 0.0f;
}

// main: Evaluate the integral I over [0,2] using qtrap, qsimp, and qromb.
// Print each result and the number of trapzd calls.
int main() {
    float a = 0.0f, b = 2.0f;

    // qtrap
    trapzdCount = 0;
    float ansQtrap = qtrap(func, a, b);
    int callsQtrap = trapzdCount;

    // qsimp
    trapzdCount = 0;
    float ansQsimp = qsimp(func, a, b);
    int callsQsimp = trapzdCount;

    // qromb: Expect convergence after exactly 5 calls.
    trapzdCount = 0;
    float ansQromb = qromb(func, a, b);
    int callsQromb = trapzdCount;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Integral I = ∫[0..2] x^4 log(x + sqrt(x^2+1)) dx\n\n";
    std::cout << "qtrap  result = " << ansQtrap  << ", trapzd calls = " << callsQtrap  << "\n";
    std::cout << "qsimp  result = " << ansQsimp  << ", trapzd calls = " << callsQsimp  << "\n";
    std::cout << "qromb  result = " << ansQromb  << ", trapzd calls = " << callsQromb  << "\n";
    
    return 0;
}
