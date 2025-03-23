// main.cpp
// assignment2_week4
// Created by Bilgin Demir on 3/12/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

#define MAXIT 10  // Maximum iterations for Newton's method
#define EPS 1e-14 // Convergence tolerance

// Gamma function approximation (needed for Gauss-Laguerre)
double gammln(double xx) {
    static double cof[6] = {76.18009172947146, -86.50532032941677, 
                            24.01409824083091, -1.231739572450155, 
                            0.1208650973866179e-2, -0.5395239384953e-5};
    double x, y, tmp, ser;
    x = y = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (int j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

// **1. Gauss-Legendre Quadrature**
void gauleg(double x1, double x2, double x[], double w[], int n) {
    int m = (n + 1) / 2;
    double xm = 0.5 * (x2 + x1);
    double xl = 0.5 * (x2 - x1);
    double z, z1, p1, p2, p3, pp;

    for (int i = 1; i <= m; i++) {
        z = cos(3.141592654 * (i - 0.25) / (n + 0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > EPS);
        x[i - 1] = xm - xl * z;
        x[n - i] = xm + xl * z;
        w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - i] = w[i - 1];
    }
}

// **2. Gauss-Chebyshev Quadrature**
void gaucheb(double x[], double w[], int n) {
    for (int i = 1; i <= n; i++) {
        x[i - 1] = cos(M_PI * (i - 0.5) / n);  // Compute Chebyshev roots
        w[i - 1] = M_PI / n;  // Weights for Chebyshev quadrature
    }
}

// **3. Gauss-Hermite Quadrature**
void gauher(double x[], double w[], int n) {
    int m = (n + 1) / 2;
    double z, z1, p1, p2, p3, pp, factorial = 1.0;
    
    // Factorial computation for weight correction
    for (int i = 1; i <= n; i++) factorial *= i;
    
    for (int i = 1; i <= m; i++) {
        if (i == 1) z = sqrt(2.0 * n + 1.0) - 1.85575 * pow(2.0 * n + 1.0, -0.16667);
        else if (i == 2) z -= 1.14 * pow(n, 0.426) / z;
        else if (i == 3) z = 1.86 * z - 0.86 * x[0];
        else if (i == 4) z = 1.91 * z - 0.91 * x[1];
        else z = 2.0 * z - x[i - 2];

        for (int its = 1; its <= MAXIT; its++) {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0 / j) * p2 - sqrt((j - 1.0) / j) * p3;
            }
            pp = sqrt(2.0 * n) * p2;
            z1 = z;
            z -= p1 / pp;
            if (fabs(z - z1) <= EPS) break;
        }

        x[i - 1] = z;
        x[n - i] = -z;
        w[i - 1] = (factorial * sqrt(M_PI)) / (n * n * pp * pp); 
        w[n - i] = w[i - 1];
    }
}

// **Main Function to Test Quadrature Computations**
int main() {
    const int n = 5;
    double x[n], w[n];

    cout << "\nGauss-Legendre Quadrature:" << endl;
    gauleg(-1, 1, x, w, n);
    for (int i = 0; i < n; i++) cout << "x[" << i << "]=" << x[i] << " w[" << i << "]=" << w[i] << endl;

    cout << "\nGauss-Chebyshev Quadrature:" << endl;
    gaucheb(x, w, n);
    for (int i = 0; i < n; i++) cout << "x[" << i << "]=" << x[i] << " w[" << i << "]=" << w[i] << endl;

    cout << "\nGauss-Hermite Quadrature:" << endl;
    gauher(x, w, n);
    for (int i = 0; i < n; i++) cout << "x[" << i << "]=" << x[i] << " w[" << i << "]=" << w[i] << endl;

    return 0;
}


