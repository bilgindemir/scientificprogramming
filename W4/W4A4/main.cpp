// main.cpp
// assignment4_week4
// Created by Bilgin Demir on 3/13/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <stdio.h>
#include <math.h>

#define EPS 1.0e-6  // Convergence tolerance for numerical integration

// Function prototypes
float qgaus(float (*func)(float), float a, float b);
float quad3d(float (*func)(float, float, float), float x1, float x2);
float H(float x);
float G(float y);
float Func(float z);

// Integration limits (replace these with desired functions)
float yy1(float x) { return -sqrt(1 - x * x); }  // Lower limit for y
float yy2(float x) { return sqrt(1 - x * x); }   // Upper limit for y
float z1(float x, float y) { return -sqrt(1 - x * x - y * y); }  // Lower limit for z
float z2(float x, float y) { return sqrt(1 - x * x - y * y); }   // Upper limit for z

// Global variables for recursive function calls
static float xsav, ysav;
static float (*nrfunc)(float, float, float);

// Example function to integrate: f(x, y, z) = x^2 + y^2 + z^2
float example_func(float x, float y, float z) {
    return x * x + y * y + z * z;
}

// Main function
int main() {
    float result;

    // Compute integral over unit sphere
    result = quad3d(example_func, -1.0, 1.0);
    
    printf("Integral result: %f\n", result);
    return 0;
}

// Top-level function for 3D integral
float quad3d(float (*func)(float, float, float), float x1, float x2) {
    nrfunc = func;  // Store function pointer globally
    return qgaus(H, x1, x2);
}

// Second integration over y
float H(float x) {
    xsav = x;  // Store x value
    return qgaus(G, yy1(x), yy2(x));  // Integrate over y
}

// Third integration over z
float G(float y) {
    ysav = y;  // Store y value
    return qgaus(Func, z1(xsav, y), z2(xsav, y));  // Integrate over z
}

// Evaluate function f(x, y, z) at current x, y, z
float Func(float z) {
    return (*nrfunc)(xsav, ysav, z);
}

// Gaussian quadrature integration (10-point Gauss-Legendre quadrature)
float qgaus(float (*func)(float), float a, float b) {
    int j;
    float xr, xm, dx, s;
    
    // 10-point Gauss-Legendre abscissas and weights
    static float x[] = {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
    static float w[] = {0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

    xm = 0.5 * (b + a);
    xr = 0.5 * (b - a);
    s = 0.0;  // Integral sum

    for (j = 1; j <= 5; j++) {
        dx = xr * x[j];
        s += w[j] * ((*func)(xm + dx) + (*func)(xm - dx));
    }
    
    return s * xr;  // Final integral value
}
