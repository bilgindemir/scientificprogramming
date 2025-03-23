// main.cpp
// assignment5_week4
// Created by Bilgin Demir on 3/16/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <stdio.h>
#include <math.h>

#define TERMS 10  // Number of terms in the series expansion

// Gamma function approximation
float gammln(float xx) {
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
                            24.01409824083091, -1.231739572450155,
                            0.1208650973866179e-2, -0.5395239384953e-5};
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

// Factorial function using Gamma function for large values
float factorial(int n) {
    static int ntop = 4;
    static float a[33] = {1.0, 1.0, 2.0, 6.0, 24.0}; // Precomputed values
    int j;
    if (n < 0) {
        printf("Negative factorial error\n");
        return -1;
    }
    if (n > 32) return exp(gammln(n + 1.0));
    while (ntop < n) { // Compute factorial iteratively
        j = ntop++;
        a[ntop] = a[j] * ntop;
    }
    return a[n];
}

// Sine series computation using Taylor series
float sin_series(float x, int terms) {
    float sum = 0.0;
    for (int j = 0; j < terms; j++) {
        float term = pow(-1, j) * pow(x, 2 * j + 1) / factorial(2 * j + 1);
        sum += term;
    }
    return sum;
}

// Euler summation to accelerate convergence of alternating series
void eulsum(float *sum, float term, int jterm, float wksp[]) {
    int j;
    static int nterm;
    float tmp, dum;

    if (jterm == 1) {
        nterm = 1;
        *sum = 0.5 * (wksp[1] = term);
    } else {
        tmp = wksp[1];
        wksp[1] = term;
        for (j = 1; j <= nterm - 1; j++) {
            dum = wksp[j + 1];
            wksp[j + 1] = 0.5 * (wksp[j] + tmp);
            tmp = dum;
        }
        wksp[nterm + 1] = 0.5 * (wksp[nterm] + tmp);
        if (fabs(wksp[nterm + 1]) <= fabs(wksp[nterm])) {
            *sum += (0.5 * wksp[++nterm]);
        } else {
            *sum += wksp[nterm + 1];
        }
    }
}

// Sine series computation with Euler summation
float sin_series_euler(float x, int terms) {
    float sum = 0.0;
    float wksp[terms + 1]; // Workspace for Euler summation
    for (int j = 0; j < terms; j++) {
        float term = pow(-1, j) * pow(x, 2 * j + 1) / factorial(2 * j + 1);
        eulsum(&sum, term, j + 1, wksp);
    }
    return sum;
}

// Main function
int main() {
    float x = 1.0;  // Input value for sin(x)
    int terms = TERMS;

    float result_standard = sin_series(x, terms);
    float result_euler = sin_series_euler(x, terms);
    float result_builtin = sin(x);

    printf("sin(%f) computed using Taylor series: %f\n", x, result_standard);
    printf("sin(%f) computed using Euler summation: %f\n", x, result_euler);
    printf("sin(%f) using built-in function: %f\n", x, result_builtin);

    return 0;
}
