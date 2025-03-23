// main.cpp
// assignment1_week4
// Created by Bilgin Demir on 3/11/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
#include <cmath>

// Define function pointers for the integrands
#define FUNC1(x) (2.0 * (x) * my_func1(a + (x) * (x)))
#define FUNC2(x) (2.0 * (x) * my_func2(b - (x) * (x)))
#define FUNC3(x) (my_func3(-log(x)) / (x))

using namespace std;

// Example function 1: Singular at lower limit (e.g., sqrt singularity)
float my_func1(float x) {
    return exp(-x * x);  // Example: Gaussian function
}

// Example function 2: Singular at upper limit
float my_func2(float x) {
    return 1.0 / sqrt(x + 1.0);  // Example: Singular function at x = -1
}

// Example function 3: Exponentially decreasing at infinity
float my_func3(float x) {
    return exp(-x);  // Example: Exponential decay
}

// Midpoint integration for inverse square-root singularity at lower limit
float midsql(float (*funk)(float), float a, float b, int n) {
    float x, tnm, sum, del, ddel;
    static float s;
    int it, j;

    float A = 0.0;
    float B = sqrt(b - a);

    if (n == 1) {
        return (s = (B - A) * FUNC1(0.5 * (A + B)));
    } else {
        for (it = 1, j = 1; j < n - 1; j++) it *= 3;
        tnm = it;
        del = (B - A) / (3.0 * tnm);
        ddel = del + del;
        x = A + 0.5 * del;
        sum = 0.0;

        for (j = 1; j <= it; j++) {
            sum += FUNC1(x);
            x += ddel;
            sum += FUNC1(x);
            x += del;
        }

        s = (s + (B - A) * sum / tnm) / 3.0;
        return s;
    }
}

// Midpoint integration for inverse square-root singularity at upper limit
float midsqu(float (*funk)(float), float a, float b, int n) {
    float x, tnm, sum, del, ddel;
    static float s;
    int it, j;

    float A = 0.0;
    float B = sqrt(b - a);

    if (n == 1) {
        return (s = (B - A) * FUNC2(0.5 * (A + B)));
    } else {
        for (it = 1, j = 1; j < n - 1; j++) it *= 3;
        tnm = it;
        del = (B - A) / (3.0 * tnm);
        ddel = del + del;
        x = A + 0.5 * del;
        sum = 0.0;

        for (j = 1; j <= it; j++) {
            sum += FUNC2(x);
            x += ddel;
            sum += FUNC2(x);
            x += del;
        }

        s = (s + (B - A) * sum / tnm) / 3.0;
        return s;
    }
}

// Midpoint integration for infinite upper limit (Exponential decay)
float midexp(float (*funk)(float), float a, float b, int n) {
    float x, tnm, sum, del, ddel;
    static float s;
    int it, j;

    float A = 0.0;
    float B = exp(-a); // Assumed b = infinity

    if (n == 1) {
        return (s = (B - A) * FUNC3(0.5 * (A + B)));
    } else {
        for (it = 1, j = 1; j < n - 1; j++) it *= 3;
        tnm = it;
        del = (B - A) / (3.0 * tnm);
        ddel = del + del;
        x = A + 0.5 * del;
        sum = 0.0;

        for (j = 1; j <= it; j++) {
            sum += FUNC3(x);
            x += ddel;
            sum += FUNC3(x);
            x += del;
        }

        s = (s + (B - A) * sum / tnm) / 3.0;
        return s;
    }
}

int main() {
    float a = 0.0, b = 1.0;
    int n = 5;  // Number of iterations

    // Test midsql (Singularity at lower bound)
    float result1 = midsql(my_func1, a, b, n);
    cout << "Integral result (midsql): " << result1 << endl;

    // Test midsqu (Singularity at upper bound)
    float result2 = midsqu(my_func2, a, b, n);
    cout << "Integral result (midsqu): " << result2 << endl;

    // Test midexp (Exponential decay)
    float result3 = midexp(my_func3, a, INFINITY, n);
    cout << "Integral result (midexp): " << result3 << endl;

    return 0;
}
