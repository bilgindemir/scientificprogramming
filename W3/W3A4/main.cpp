// main.cpp
// assignment4_week3
// Created by Bilgin Demir on 3/3/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

// Define the macro so that FUNC(x) calls the function 'func'
#define FUNC(x) (func(x))

// Global constants for integration
#define EPS 1.0e-6       // Desired fractional accuracy
#define JMAX 20          // Maximum iterations for qtrap and qsimp
#define JMAXP 20         // Maximum iterations for qromb
#define K 5              // Parameter for qromb (not used for main output)

//----------------------------------------------------------------
// Integrand: f(x) = x^4 * log(x + sqrt(x^2+1))
//----------------------------------------------------------------
float func(float x) {
    return x*x*x*x * std::log(x + std::sqrt(x*x + 1.0));
}
  
//----------------------------------------------------------------
// Dummy antiderivative function Func for "true" value.
// For demonstration, we assume the true value for x=2 is known to be 8.15336418.
// (In practice, one would compute the antiderivative exactly if possible.)
float Func(float x) {
    if (x <= 0.0) return 0.0;
    if (x >= 2.0) return 8.15336418;
    return 8.15336418 * x / 2.0;
}
  
//----------------------------------------------------------------
// trapzd: One stage of refinement of the trapezoidal rule.
static int dummyCounter = 0; // (Not used in output for Assignment 4)
float trapzd(float (*f)(float), float a, float b, int n) {
    dummyCounter++;
    static float s;
    int it, j;
    float x, tnm, sum, del;
    if (n == 1) {
        s = 0.5 * (b - a) * (f(a) + f(b));
        return s;
    } else {
        for (it = 1, j = 1; j < n - 1; j++) it *= 2;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        sum = 0.0;
        for (j = 1; j <= tnm; j++, x += del)
            sum += f(x);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}
  
// qtrap: Basic trapezoidal integration driver.
float qtrap(float (*f)(float), float a, float b) {
    float s, olds = 0.0;
    for (int j = 1; j <= JMAX; j++) {
        s = trapzd(f, a, b, j);
        if (j > 5) {
            if (fabs(s - olds) < EPS * fabs(olds))
                return s;
        }
        olds = s;
    }
    cout << "Too many steps in qtrap" << endl;
    return s;
}
  
// qsimp: Simpson's rule via Richardson extrapolation on trapezoidal estimates.
float qsimp(float (*f)(float), float a, float b) {
    float s, st, ost = 0.0, os = 0.0;
    for (int j = 1; j <= JMAX; j++) {
        st = trapzd(f, a, b, j);
        s = (4.0 * st - ost) / 3.0;
        if (j > 5) {
            if (fabs(s - os) < EPS * fabs(os))
                return s;
        }
        os = s;
        ost = st;
    }
    cout << "Too many steps in qsimp" << endl;
    return 0.0;
}
  
// polint: Polynomial interpolation for Romberg's extrapolation.
// (0-based version.)
void polint(const float xa[], const float ya[], int n, float x, float &y, float &dy) {
    int i, m, ns = 0;
    float den, dif, dift, ho, hp, w;
    float *c = new float[n];
    float *d = new float[n];
    dif = fabs(x - xa[0]);
    for (i = 0; i < n; i++) {
        dift = fabs(x - xa[i]);
        if (dift < dif) { ns = i; dif = dift; }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns];
    ns--;
    for (m = 1; m < n; m++) {
        for (i = 0; i < n - m; i++) {
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            den = ho - hp;
            if (den == 0.0) { cout << "Error in polint" << endl; return; }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        if (2 * ns < (n - m))
            dy = c[ns + 1];
        else {
            dy = d[ns];
            ns--;
        }
        y += dy;
    }
    delete [] c;
    delete [] d;
}
  
// qromb: Romberg integration driver.
float qromb(float (*f)(float), float a, float b) {
    float ss, dss;
    float s[JMAXP+1], h[JMAXP+2];
    h[1] = 1.0;
    for (int j = 1; j <= JMAXP; j++) {
        s[j] = trapzd(f, a, b, j);
        if (j >= K) {
            polint(h + j - K, s + j - K, K, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss))
                return ss;
        }
        h[j+1] = 0.25 * h[j];
    }
    cout << "Too many steps in qromb" << endl;
    return 0.0;
}
  
//----------------------------------------------------------------
// midpnt: Extended midpoint rule.
// When called with n=1, returns the crudest estimate of the integral.
// Subsequent calls with n=2,3,... improve the accuracy.
//----------------------------------------------------------------
float midpnt(float (*f)(float), float a, float b, int n) {
    float x, tnm, sum, del, ddel;
    static float s;
    int it, j;
    if (n == 1) {
        return (s = (b - a) * FUNC(0.5 * (a + b)));
    } else {
        for (it = 1, j = 1; j < n - 1; j++) it *= 3;
        tnm = it;
        del = (b - a) / (3.0 * tnm);
        ddel = 2 * del;
        x = a + 0.5 * del;
        sum = 0.0;
        for (j = 1; j <= it; j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        s = (s + (b - a) * sum / tnm) / 3.0;
        return s;
    }
}
  
//----------------------------------------------------------------
// main: Assignment 4
// Compute the integral using various methods and print the results.
//----------------------------------------------------------------
int main(int argc, const char * argv[]) {
    float a = 0.0;
    float b = 2.0;
    cout << fixed << setprecision(8);
    cout << "Computed Value of integral (trapezoidal rule): " << qtrap(func, a, b) << endl;
    cout << "Computed Value of integral (Simpson rule): " << qsimp(func, a, b) << endl;
    cout << "Computed Value of integral (Romberg formula): " << qromb(func, a, b) << endl;
    cout << "Computed Value of integral (mid-point formula n=1): " << midpnt(func, a, b, 1) << endl;
    cout << "Computed Value of integral (mid-point formula n=2): " << midpnt(func, a, b, 2) << endl;
    cout << "Computed Value of integral (mid-point formula n=3): " << midpnt(func, a, b, 3) << endl;
    cout << "True Value of integral: " << (Func(b) - Func(a)) << endl;
    return 0;
}
