// main.cpp
// assignment2_week3
// Created by Bilgin Demir on 3/3/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>
#include <cstdlib>

//-------------------------------
// Set desired parameters:
#define EPS 1.0e-6   // Desired fractional accuracy
#define JMAX 20      // Maximum allowed iterations

// Macro to call the integrand (for clarity)
#define FUNC(x) (func(x))

//-------------------------------
// 1) Define the integrand: f(x) = sin(x)
//    (You can change this function as needed)
float func(float x) {
    return std::sin(x);
}

//-------------------------------
// 2) trapzd: One stage of refinement of the trapezoidal rule.
//    When n==1, it returns the crude estimate. For n>1, it adds additional points.
float trapzd(float (*funk)(float), float a, float b, int n) {
    static float s; // holds the previous trapezoidal estimate
    int it, j;
    float x, tnm, sum, del;

    if (n == 1) {
        s = 0.5f * (b - a) * (FUNC(a) + FUNC(b));
        return s;
    } else {
        // Determine the number of new points to add: it = 2^(n-2)
        for (it = 1, j = 1; j < n - 1; j++) {
            it <<= 1;  // equivalent to it = it * 2
        }
        tnm = it;
        del = (b - a) / tnm;  // spacing between new points
        x = a + 0.5f * del;
        sum = 0.0f;
        for (j = 1; j <= it; j++, x += del) {
            sum += FUNC(x);
        }
        s = 0.5f * (s + (b - a) * sum / tnm);
        return s;
    }
}

// 3) qsimp: Simpson’s rule using extrapolation from trapzd.
//    It uses the formula s = (4*st - ost)/3, comparing successive estimates.
float qsimp(float (*func)(float), float a, float b) {
    int j;
    float s, st, ost = 0.0f, os = 0.0f;
    for (j = 1; j <= JMAX; j++) {
        st = trapzd(func, a, b, j);
        s = (4.0f * st - ost) / 3.0f;  // Simpson's rule extrapolation
        if (j > 5) {  // Avoid spurious early convergence.
            if ( (std::fabs(s - os) < EPS * std::fabs(os)) || 
                 (s == 0.0f && os == 0.0f) )
                return s;
        }
        os = s;
        ost = st;
    }
    std::cout << "Too many steps in routine qsimp" << std::endl;
    return 0.0f; // Should not normally get here.
}

//-------------------------------
// 4) main: Test Simpson’s rule integration
int main() {
    float a = 0.0f;          // Lower limit of integration
    float b = M_PI;          // Upper limit of integration (integrate sin(x) from 0 to pi)
    
    float result = qsimp(func, a, b);
    
    std::cout << "Simpson's rule integration of sin(x) from " 
              << a << " to " << b << " = " << result << std::endl;
    
    // Exact integral of sin(x) over [0, pi] is 2
    float exact = 2.0f;
    std::cout << "Exact value: " << exact << std::endl;
    std::cout << "Absolute error: " << std::fabs(result - exact) << std::endl;
    
    return 0;
}
