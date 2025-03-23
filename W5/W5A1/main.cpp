// main.cpp
// assignment1_week5
// Created by Bilgin Demir on 3/23/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.

#include <iostream>
using namespace std;

// Evaluate polynomial and derivatives
void ddpoly(float c[], int nc, float x, float pd[], int nd) {
    int nnd, j, i;
    float cnst = 1.0;
    pd[0] = c[nc];
    for (j = 1; j <= nd; j++) pd[j] = 0.0;

    for (i = nc - 1; i >= 0; i--) {
        nnd = (nd < (nc - i) ? nd : nc - i);
        for (j = nnd; j >= 1; j--)
            pd[j] = pd[j] * x + pd[j - 1];
        pd[0] = pd[0] * x + c[i];
    }

    for (i = 2; i <= nd; i++) {
        cnst *= i;
        pd[i] *= cnst;
    }
}

// Polynomial division
void poldiv(float u[], int n, float v[], int nv, float q[], float r[]) {
    int k, j;
    for (j = 0; j <= n; j++) {
        r[j] = u[j];
        q[j] = 0.0;
    }
    for (k = n - nv; k >= 0; k--) {
        q[k] = r[nv + k] / v[nv];
        for (j = nv + k - 1; j >= k; j--)
            r[j] -= q[k] * v[j - k];
    }
    for (j = nv; j <= n; j++) r[j] = 0.0;
}

int main() {
    const int MAX = 100;
    float c[MAX], pd[MAX];
    int degree, nd;
    float x;

    // Part 1: Evaluate polynomial and derivatives
    cout << "Enter degree N of polynomial P(x): ";
    cin >> degree;

    cout << "Enter " << degree + 1 << " coefficients (c[0] to c[" << degree << "]):\n";
    for (int i = 0; i <= degree; i++) cin >> c[i];

    cout << "Enter value of x: ";
    cin >> x;

    cout << "Enter number of derivatives to compute: ";
    cin >> nd;

    ddpoly(c, degree, x, pd, nd);

    cout << "\nEvaluation at x = " << x << ":\n";
    cout << "P(x) = " << pd[0] << endl;
    for (int i = 1; i <= nd; i++)
        cout << "P^(" << i << ")(x) = " << pd[i] << endl;

    // Part 2: Polynomial division
    int deg_u, deg_v;
    float u[MAX], v[MAX], q[MAX], r[MAX];

    cout << "\n--- Polynomial Division ---\n";
    cout << "Enter degree of numerator polynomial u(x): ";
    cin >> deg_u;
    cout << "Enter " << deg_u + 1 << " coefficients of u(x):\n";
    for (int i = 0; i <= deg_u; i++) cin >> u[i];

    cout << "Enter degree of denominator polynomial v(x): ";
    cin >> deg_v;
    cout << "Enter " << deg_v + 1 << " coefficients of v(x):\n";
    for (int i = 0; i <= deg_v; i++) cin >> v[i];

    poldiv(u, deg_u, v, deg_v, q, r);

    cout << "\nQuotient polynomial q(x):\n";
    for (int i = 0; i <= deg_u; i++) cout << "q[" << i << "] = " << q[i] << endl;

    cout << "\nRemainder polynomial r(x):\n";
    for (int i = 0; i <= deg_u; i++) cout << "r[" << i << "] = " << r[i] << endl;

    return 0;
}