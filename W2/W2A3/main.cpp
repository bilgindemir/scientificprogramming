//
// main.cpp
// assignment3_week2
//
// Created by Bilgin Demir on 2/27/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.
// 

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
using namespace std;

// ---------------------------------------------------------------------------
// We store the final interpolation result in a global variable, yvalue.
// This matches the style of the Numerical Recipes splint routine.
static float yvalue;  

// ---------------------------------------------------------------------------
// Function to be interpolated.
// Example: f(x) = 3x^2 + (1 / π^4) * log((π - x)^2) + 1
float fun(float x) {
    return 3.0f * x * x 
         + (1.0f / (float)pow(M_PI, 4)) * log((M_PI - x)*(M_PI - x)) 
         + 1.0f;
}

// ---------------------------------------------------------------------------
// A simple selection sort to keep x-array in ascending order (1-indexed).
void selectionSort(int start, int end, float arr[]) {
    for (int i = start; i < end; i++) {
        int minIndex = i;
        for (int j = i + 1; j <= end; j++) {
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }
        float temp = arr[i];
        arr[i] = arr[minIndex];
        arr[minIndex] = temp;
    }
}

// ---------------------------------------------------------------------------
// Numerical Recipes: spline(x[], y[], n, yp1, ypn, y2[])
//  - x[1..n], y[1..n]: Tabulated function (x[i], y[i])
//  - n: Number of points
//  - yp1, ypn: First derivatives at the boundaries (or large number for "natural" spline)
//  - y2[]: Output array of second derivatives
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]) {
    int i, k;
    float p, qn, sig, un;
    // Temporary array for the tridiagonal system
    float u[n+1];

    // Boundary condition at x[1]
    if (yp1 > 0.99e+30f) {
        // Natural spline
        y2[1] = 0.0f;
        u[1]  = 0.0f;
    } else {
        // Specified first derivative
        y2[1] = -0.5f;
        u[1]  = (3.0f/(x[2] - x[1])) * ((y[2] - y[1])/(x[2] - x[1]) - yp1);
    }

    // Decomposition loop of the tridiagonal algorithm
    for (i = 2; i <= n - 1; i++) {
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        p   = sig * y2[i-1] + 2.0f;
        y2[i] = (sig - 1.0f) / p;
        u[i]  = ( (y[i+1] - y[i]) / (x[i+1] - x[i])
                - (y[i] - y[i-1]) / (x[i] - x[i-1]) );
        u[i]  = (6.0f * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
    }

    // Boundary condition at x[n]
    if (ypn > 0.99e+30f) {
        // Natural spline
        qn = 0.0f;
        un = 0.0f;
    } else {
        // Specified first derivative
        qn = 0.5f;
        un = (3.0f/(x[n] - x[n-1])) * (ypn - (y[n] - y[n-1])/(x[n] - x[n-1]));
    }

    y2[n] = (un - qn * u[n-1]) / (qn * y2[n-1] + 1.0f);

    // Backsubstitution loop
    for (k = n - 1; k >= 1; k--) {
        y2[k] = y2[k] * y2[k+1] + u[k];
    }
}

// ---------------------------------------------------------------------------
// Numerical Recipes: splint(xa[], ya[], y2a[], n, x)
//  - xa[1..n], ya[1..n]: Tabulated function (already used in spline)
//  - y2a[1..n]: Second derivatives from spline()
//  - n: Number of points
//  - x: The value at which we want to interpolate
// On output, the global yvalue is set to the spline interpolation result.
void splint(float xa[], float ya[], float y2a[], int n, float x) {
    int klo = 1;
    int khi = n;

    // Bisection search to find the correct interval [klo, khi]
    while (khi - klo > 1) {
        int k = (khi + klo) >> 1; // midpoint
        if (xa[k] > x) khi = k;
        else           klo = k;
    }

    float h = xa[khi] - xa[klo];
    if (h == 0.0f) {
        cerr << "Bad xa input to routine splint" << endl;
        return;
    }

    float a = (xa[khi] - x) / h;
    float b = (x - xa[klo]) / h;
    // Cubic spline polynomial evaluation
    yvalue = a * ya[klo] 
           + b * ya[khi] 
           + ( (a*a*a - a) * y2a[klo]
             + (b*b*b - b) * y2a[khi] ) * (h*h) / 6.0f;
}

// ---------------------------------------------------------------------------
// Main program:
//  1. Generate n data points x[i], y[i] = fun(x[i]) in ascending order.
//  2. Compute the second derivatives array y2[] using spline() with "natural" boundary conditions.
//  3. Pick a point x0 inside [x[1], x[n]] for interpolation.
//  4. Call splint() to get the interpolated value yvalue.
//  5. Compare to the "true" function value fun(x0).
// ---------------------------------------------------------------------------
int main() {
    srand(static_cast<unsigned>(time(nullptr)));

    // Number of points
    const int n = 10;
    // 1-based arrays
    float x[n+1], y[n+1], y2[n+1];

    // Generate random x-values in [0, 0.5 / π] (example)
    for (int i = 1; i <= n; i++) {
        x[i] = (float)rand() / RAND_MAX;    // in [0,1]
        x[i] *= 0.5f * (1.0f / M_PI);       // scale to [0, 0.5/π]
    }
    // Sort x
    selectionSort(1, n, x);

    // Compute y = fun(x)
    for (int i = 1; i <= n; i++) {
        y[i] = fun(x[i]);
    }

    // Print tabulated points
    cout << "Tabulated data points (x[i], y[i]):\n";
    for (int i = 1; i <= n; i++) {
        cout << "  x[" << i << "] = " << x[i] 
             << ", y[" << i << "] = " << y[i] << "\n";
    }

    // Boundary conditions: Natural spline => pass large value (e.g. 1.0e30)
    float yp1 = 1.0e30f;  // derivative at x[1]
    float ypn = 1.0e30f;  // derivative at x[n]

    // Compute second derivatives
    spline(x, y, n, yp1, ypn, y2);

    // Choose x0 for interpolation within [x[1], x[n]]
    float xMin = x[1];
    float xMax = x[n];
    float x0 = xMin + ((float)rand() / RAND_MAX) * (xMax - xMin);

    cout << "\nInterpolation point x0 = " << x0 << " (within [" 
         << xMin << ", " << xMax << "])\n";

    // Perform cubic spline interpolation
    splint(x, y, y2, n, x0);

    // Compare to the "true" function value
    float trueVal = fun(x0);

    cout << "\nResults:\n";
    cout << "  x0           = " << x0 << "\n"
         << "  True f(x0)   = " << trueVal << "\n"
         << "  Spline value = " << yvalue << "\n"
         << "  Error        = " << (yvalue - trueVal) << "\n\n";

    return 0;
}
