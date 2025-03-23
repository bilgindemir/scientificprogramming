// main.cpp
// assignment3_week5
// Created by Bilgin Demir on 3/23/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

// Constants
#define NTAB 10
#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define SAFE 2.0

// Sample function to differentiate (you can change this)
float sample_func(float x) {
    return x * x;  // Derivative should be 2x, for testing
}

// Ridders' method for numerical differentiation
float dfridr(float (*func)(float), float x, float h, float *err) {
    int i, j;
    float errt, fac, hh, ans = 0;
    float a[NTAB+1][NTAB+1];

    if (h == 0.0) {
        cout << "h must be nonzero in dfridr." << endl;
        exit(1);
    }

    hh = h;
    a[1][1] = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
    *err = BIG;

    for (i = 2; i <= NTAB; i++) {
        hh /= CON;
        a[1][i] = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
        fac = CON2;

        for (j = 2; j <= i; j++) {
            a[j][i] = (a[j-1][i]*fac - a[j-1][i-1]) / (fac - 1.0);
            fac = CON2 * fac;
            errt = fmax(fabs(a[j][i] - a[j-1][i]), fabs(a[j][i] - a[j-1][i-1]));
            if (errt <= *err) {
                *err = errt;
                ans = a[j][i];
            }
        }

        if (fabs(a[i][i] - a[i-1][i-1]) >= SAFE * (*err)) break;
    }

    return ans;
}

int main() {
    float x, h, err;

    cout << "Numerical Derivative using Ridders' Method\n";
    cout << "------------------------------------------\n";
    cout << "We will compute the derivative of f(x) = x^2\n";

    cout << "Enter the value of x: ";
    cin >> x;

    cout << "Enter initial step size h (e.g., 0.1): ";
    cin >> h;

    float result = dfridr(sample_func, x, h, &err);

    cout << "\nEstimated derivative at x = " << x << " is: " << result << endl;
    cout << "Estimated error: " << err << endl;

    return 0;
}
