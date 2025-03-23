// main.cpp
// assignment2_week5
// Created by Bilgin Demir on 3/23/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
using namespace std;

// Evaluate rational function R(x) = P(x)/Q(x)
// P(x) = cof[0] + cof[1]x + ... + cof[mm]x^mm
// Q(x) = 1 + cof[mm+1]x + ... + cof[mm+kk]x^kk
float ratval(float x, float cof[], int mm, int kk) {
    int j;
    float sumd, sumn;

    // Evaluate numerator P(x)
    sumn = cof[mm];
    for (j = mm - 1; j >= 0; j--) {
        sumn = sumn * x + cof[j];
    }

    // Evaluate denominator Q(x), excluding the leading 1
    sumd = 0.0;
    for (j = mm + kk; j >= mm + 1; j--) {
        sumd = (sumd + cof[j]) * x;
    }

    return sumn / (1.0 + sumd);
}

int main() {
    const int MAX = 100;
    float cof[MAX];
    int mm, kk;
    float x;

    cout << "Enter degree of numerator polynomial P(x): ";
    cin >> mm;

    cout << "Enter degree of denominator polynomial Q(x): ";
    cin >> kk;

    cout << "Enter " << mm + 1 << " coefficients for numerator (p[0] to p[" << mm << "]):\n";
    for (int i = 0; i <= mm; i++) {
        cin >> cof[i];
    }

    cout << "Enter " << kk << " coefficients for denominator (q[1] to q[" << kk << "]), assuming q[0] = 1:\n";
    for (int i = 1; i <= kk; i++) {
        cin >> cof[mm + i];
    }

    cout << "Enter value of x to evaluate R(x): ";
    cin >> x;

    float result = ratval(x, cof, mm, kk);

    cout << "\nR(x) = P(x)/Q(x) evaluated at x = " << x << " is: " << result << endl;

    return 0;
}
