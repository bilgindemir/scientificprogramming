// main.cpp
// assignment2_week4
// Created by Bilgin Demir on 3/11/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>

using namespace std;

// Define the function to be integrated
float my_function(float x) {
    return exp(-x * x);  // Example: Gaussian function
}

// Gauss-Legendre Quadrature Integration (N=10)
float qgaus(float (*func)(float), float a, float b) {
    int j;
    float xr, xm, dx, s;

    // Gauss-Legendre 10-point abscissas (x) and weights (w)
    static float x[] = {0.0, 0.1488743389, 0.4333953941, 
                        0.6794095682, 0.8650633666, 0.9739065285};
    static float w[] = {0.0, 0.2955242247, 0.2692667193, 
                        0.2190863625, 0.1494513491, 0.0666713443};

    // Midpoint and half-range
    xm = 0.5 * (b + a);
    xr = 0.5 * (b - a);
    s = 0.0;

    // Compute the weighted sum
    for (j = 1; j <= 5; j++) {
        dx = xr * x[j];
        s += w[j] * ((*func)(xm + dx) + (*func)(xm - dx));
    }

    return s * xr; // Scale to the integration range
}

int main() {
    float a = -1.0, b = 1.0;  // Integration limits
    float result = qgaus(my_function, a, b);
    
    cout << "Integral result (Gauss-Legendre Quadrature): " << result << endl;
    return 0;
}
