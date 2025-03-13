// main.cpp
// assignment1_week3
// Created by Bilgin Demir on 3/3/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>
#include <cstdlib>

// These can be adjusted as needed:
#define EPS 1.0e-6   // Desired fractional accuracy
#define JMAX 20      // Maximum allowed iterations (2^(JMAX-1) subintervals)

// Macro so that trapzd matches the style from Numerical Recipes
#define FUNC(x) (func(x))

// 1) Define the function to be integrated.
//    For demonstration, I use f(x) = sin(x). Replace with any function you want.
float func(float x) {
    return std::sin(x);
}

// 2) trapzd: One stage of refinement of the trapezoidal rule
//    - func is the integrand
//    - a, b are the integration limits
//    - n is which refinement step we are on
//    - The static variable s stores the previous trapezoid approximation.
float trapzd(float (*funk)(float), float a, float b, int n) {
    static float s;
    int it, j;
    float x, tnm, sum, del;

    if (n == 1) {
        // Crudest estimate (n=1): area of one trapezoid with corners at a,b
        s = 0.5f * (b - a) * (FUNC(a) + FUNC(b));
        return s;
    }
    else {
        // At each subsequent call, we add 2^(n-2) equally spaced interior points
        for (it = 1, j = 1; j < n - 1; j++) {
            it <<= 1;  // 2^(n-2)
        }
        tnm = it;
        del = (b - a) / tnm;   // Spacing between new points
        x = a + 0.5f * del;

        sum = 0.0f;
        for (j = 1; j <= it; j++, x += del) {
            sum += FUNC(x);
        }
        // Refine the previous trapezoid value s
        s = 0.5f * (s + (b - a) * sum / tnm);
        return s;
    }
}

// 3) qtrap: High-level driver using repeated calls to trapzd.
//    It returns the integral of the function func from a to b
//    to within a fractional accuracy EPS, or until JMAX refinements are used.
float qtrap(float (*func)(float), float a, float b) {
    float s, olds = 0.0f;
    // Perform successive trapezoid refinements until convergence or JMAX steps
    for (int j = 1; j <= JMAX; j++) {
        s = trapzd(func, a, b, j);

        // After a few initial steps, check for convergence
        if (j > 5) {
            if ((std::fabs(s - olds) < EPS * std::fabs(olds)) || 
                (s == 0.0f && olds == 0.0f)) {
                return s;
            }
        }
        olds = s;
    }
    std::cerr << "Too many steps in routine qtrap" << std::endl;
    return s; // Fallback (should rarely happen if JMAX is large enough)
}

// 4) main: Demonstration of trapezoidal integration
int main() {
    float a = 0.0f;  // lower limit
    float b = M_PI;  // upper limit (example: integrate sin(x) from 0 to π)

    float result = qtrap(func, a, b);
    std::cout << "Integral of sin(x) from " << a << " to " << b << " = " 
              << result << std::endl;

    // Compare to the exact answer = 2 for sin(x) from 0..π
    float exact = 2.0f;
    std::cout << "Exact answer: " << exact << std::endl;
    std::cout << "Error: " << (result - exact) << std::endl;

    return 0;
}
