// main.cpp
// assignment5_week3
// Created by Bilgin Demir on 3/3/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

// Global constants for the routines
#define EPS 1.0e-6   // Desired fractional accuracy
#define JMAX 10      // Maximum number of steps (kept modest due to step tripling)
#define K 4          // Number of points for polynomial extrapolation in qromo

// Prototype for polynomial interpolation (polint)
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);

// Example integrand: f(x) = x^2
float func(float x) {
    return x * x;
}

// Antiderivative of f(x)=x^2: F(x)=x^3/3 (for checking the true value)
float F(float x) {
    return x * x * x / 3.0;
}

// midpnt_re: A reentrant version of the midpoint rule using step tripling.
// It computes the nth stage of refinement of the integral of func from a to b.
float midpnt_re(float (*func)(float), float a, float b, int n) {
    float s = (b - a) * func(0.5 * (a + b));  // crude estimate for n == 1
    if (n == 1)
        return s;

    // Refine the estimate iteratively
    for (int k = 2; k <= n; k++) {
        int it = 1;
        for (int j = 1; j < k - 1; j++)
            it *= 3;  // Tripling the number of groups each stage
        float tnm = (float)it;
        float del = (b - a) / (3.0 * tnm);
        float ddel = 2.0 * del;
        float x = a + 0.5 * del;
        float sum = 0.0;
        for (int j = 1; j <= it; j++) {
            sum += func(x);
            x += ddel;
            sum += func(x);
            x += del;
        }
        s = (s + (b - a) * sum / tnm) / 3.0;
    }
    return s;
}

// Modified qtrap: Uses midpnt_re instead of trapzd.
// This routine refines the integration until successive estimates differ by less than EPS.
float qtrap(float (*func)(float), float a, float b) {
    float s, olds = 0.0;
    for (int j = 1; j <= JMAX; j++) {
        s = midpnt_re(func, a, b, j);
        if (j > 5) {  // Avoid spurious early convergence
            if (fabs(s - olds) < EPS * fabs(olds) || (s == 0.0 && olds == 0.0))
                return s;
        }
        olds = s;
    }
    cout << "Too many steps in routine qtrap" << endl;
    return 0.0;  // Should never get here.
}

// Modified qsimp: Simpson’s rule extrapolation using midpnt_re.
// The extrapolation now uses s = (9.0 * st - ost)/8.0 to reflect error reduction by 1/9.
float qsimp(float (*func)(float), float a, float b) {
    float s, st, ost = 0.0, os = 0.0;
    for (int j = 1; j <= JMAX; j++) {
        st = midpnt_re(func, a, b, j);
        s = (9.0 * st - ost) / 8.0;
        if (j > 5) {  // Avoid spurious early convergence
            if (fabs(s - os) < EPS * fabs(os) || (s == 0.0 && os == 0.0))
                return s;
        }
        os = s;
        ost = st;
    }
    cout << "Too many steps in routine qsimp" << endl;
    return 0.0;  // Should never get here.
}

// Generalized Romberg integration (qromo):
// This routine uses any integration routine "choose" (here, midpnt_re) that refines by tripling the steps.
// It uses polynomial extrapolation (via polint) to accelerate convergence.
float qromo(float (*func)(float), float a, float b,
            float (*choose)(float (*)(float), float, float, int)) {
    // Use 0-based indexing for arrays
    const int JM = JMAX;  // We'll use JMAX iterations
    float s[JM], h[JM+1];
    h[0] = 1.0;
    float ss, dss;
    for (int j = 0; j < JM; j++) {
        s[j] = choose(func, a, b, j + 1);
        if (j >= (K - 1)) {  // Once we have K values, attempt extrapolation
            // Use polynomial interpolation on the last K estimates.
            polint(&h[j - K + 1], &s[j - K + 1], K, 0.0, &ss, &dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j + 1] = h[j] / 9.0;  // Adjust step size by factor 1/9 per refinement.
    }
    cout << "Too many steps in routine qromo" << endl;
    return 0.0;  // Should never get here.
}

// polint: Polynomial interpolation based on Numerical Recipes.
// Given arrays xa and ya of length n, and a value x, it returns the interpolated value y and error dy.
void polint(float xa[], float ya[], int n, float x, float *y, float *dy) {
    int i, m, ns = 0;
    float diff = fabs(x - xa[0]);
    float *c = new float[n];
    float *d = new float[n];
    
    for (i = 0; i < n; i++) {
        float dift = fabs(x - xa[i]);
        if (dift < diff) {
            ns = i;
            diff = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    
    *y = ya[ns];
    ns = ns - 1;
    for (m = 1; m < n; m++) {
        for (i = 0; i < n - m; i++) {
            float ho = xa[i] - x;
            float hp = xa[i + m] - x;
            float w = c[i + 1] - d[i];
            float den = ho - hp;
            if (den == 0.0) {
                cerr << "Error in polint: division by zero" << endl;
                exit(1);
            }
            w /= den;
            c[i] = ho * w;
            d[i] = hp * w;
        }
        if (2 * (ns + 1) < (n - m))
            *dy = c[ns + 1];
        else {
            *dy = d[ns];
            ns--;
        }
        *y += *dy;
    }
    
    delete[] c;
    delete[] d;
}

int main() {
    float a = 0.0;
    float b = 2.0;

    cout << "Integration of f(x)=x^2 from " << a << " to " << b << "\n";
    cout << "True value (F(b)-F(a)) = " << F(b) - F(a) << "\n\n";

    cout << "Computed value using qtrap (with midpnt_re): " << qtrap(func, a, b) << "\n";
    cout << "Computed value using qsimp (with modified extrapolation): " << qsimp(func, a, b) << "\n";
    cout << "Computed value using qromo (with midpnt_re as choose): " << qromo(func, a, b, midpnt_re) << "\n";

    // Demonstrate individual stages of midpnt_re
    cout << "\nMidpoint rule results (midpnt_re):\n";
    for (int j = 1; j <= 4; j++) {
        cout << "midpnt_re(func, a, b, " << j << ") = " << midpnt_re(func, a, b, j) << "\n";
    }
    
    return 0;
}
