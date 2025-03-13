//
// main.cpp
// assignment5_week2
//
// Created by Bilgin Demir on 2/27/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.
// 

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

// -----------------------------------------------------------------------------
// We need a global or static 2D array `c[1..4][1..4]` for the coefficients,
// and global variables ansy, ansy1, ansy2 for the interpolated function and its
// partial derivatives, as per Numerical Recipes convention.
static float c[5][5];         // We'll store the bicubic coefficients here.
static float ansy, ansy1, ansy2;  // Interpolated values: f(x,y), df/dx, df/dy.

// -----------------------------------------------------------------------------
// Function f(x, y) = sin(x) * cos(y)
float fxy(float x, float y) {
    return sinf(x) * cosf(y);
}

// Partial derivatives:
// df/dx = cos(x)*cos(y)
// df/dy = -sin(x)*sin(y)
// d^2 f/(dx dy) = -cos(x)*sin(y)

// -----------------------------------------------------------------------------
// bcucof: Compute the bicubic coefficients c[1..4][1..4]
// from the function values and partial derivatives at the four corners.
void bcucof(float y[], float y1[], float y2[], float y12[], float d1, float d2) {
    // The stored weight matrix from Numerical Recipes
    static int wt[16][16] = {
        { 1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 },
        { 0,0,0,0, 0,0,0,0, 1,0,0,0, 0,0,0,0 },
        { -3,0,0,3, 0,0,0,0, -2,0,0,-1, 0,0,0,0 },
        { 2,0,0,-2, 0,0,0,0, 1,0,0,1, 0,0,0,0 },
        { 0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0 },
        { 0,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,0,0 },
        { 0,0,0,0, -3,0,0,3, 0,0,0,0, -2,0,0,-1 },
        { 0,0,0,0, 2,0,0,-2, 0,0,0,0, 1,0,0,1 },
        { -3,3,0,0, -2,-1,0,0, 0,0,0,0, 0,0,0,0 },
        { 0,0,0,0, 0,0,0,0, -3,3,0,0, -2,-1,0,0 },
        { 9,-9,9,-9, 6,3,-3,-6, 6,-6,-3,3, 4,2,1,2 },
        { -6,6,-6,6, -4,-2,2,4, -3,3,3,-3, -2,-1,-1,-2 },
        { 2,-2,0,0, 1,1,0,0, 0,0,0,0, 0,0,0,0 },
        { 0,0,0,0, 0,0,0,0, 2,-2,0,0, 1,1,0,0 },
        { -6,6,-6,6, -3,-3,3,3, -4,4,2,-2, -2,-2,-1,-1 },
        { 4,-4,4,-4, 2,2,-2,-2, 2,-2,-2,2, 1,1,1,1 }
    };

    float d1d2 = d1 * d2;
    float x[16], cl[16];

    // Pack the function values and derivatives into x[]
    // The corners are numbered 1..4 in counterclockwise order:
    //   y[1], y[2], y[3], y[4]  => function values
    //   y1[i] => partial wrt x
    //   y2[i] => partial wrt y
    //   y12[i] => cross partial wrt x,y
    x[0]  = y[1];
    x[1]  = y[2];
    x[2]  = y[3];
    x[3]  = y[4];
    x[4]  = y1[1]*d1;
    x[5]  = y1[2]*d1;
    x[6]  = y1[3]*d1;
    x[7]  = y1[4]*d1;
    x[8]  = y2[1]*d2;
    x[9]  = y2[2]*d2;
    x[10] = y2[3]*d2;
    x[11] = y2[4]*d2;
    x[12] = y12[1]*d1d2;
    x[13] = y12[2]*d1d2;
    x[14] = y12[3]*d1d2;
    x[15] = y12[4]*d1d2;

    // Multiply by the stored weight matrix
    for (int i = 0; i < 16; i++) {
        float xx = 0.0f;
        for (int k = 0; k < 16; k++) {
            xx += wt[i][k] * x[k];
        }
        cl[i] = xx;
    }

    // Unpack cl[] into the global c[1..4][1..4]
    int l = 0;
    for (int i = 1; i <= 4; i++) {
        for (int j = 1; j <= 4; j++) {
            c[i][j] = cl[l++];
        }
    }
}

// -----------------------------------------------------------------------------
// bcuint: Bicubic interpolation within a single grid cell
// y, y1, y2, y12 are function values & partial derivatives at corners
// x1l,x1u,x2l,x2u define the cell, x1,x2 is the point to interpolate.
// ansy, ansy1, ansy2 (global) get the function value and partial derivatives.
// -----------------------------------------------------------------------------
void bcuint(float y[], float y1[], float y2[], float y12[],
            float x1l, float x1u, float x2l, float x2u, float x1, float x2)
{
    // 1) Compute cell sizes
    float d1 = x1u - x1l;
    float d2 = x2u - x2l;

    if (d1 == 0.0f || d2 == 0.0f) {
        cerr << "Bad input in bcuint: zero cell size." << endl;
        ansy = ansy1 = ansy2 = 0.0f;
        return;
    }

    // 2) Compute coefficients
    bcucof(y, y1, y2, y12, d1, d2);

    // 3) Evaluate the bicubic polynomial
    float t = (x1 - x1l) / d1; // "Normalized" coordinate in x
    float u = (x2 - x2l) / d2; // "Normalized" coordinate in y

    ansy = ansy1 = ansy2 = 0.0f;

    // Evaluate from highest power to lowest (3rd-degree polynomial in t, u)
    for (int i = 4; i >= 1; i--) {
        ansy  = t*ansy  + ((c[i][4]*u + c[i][3])*u + c[i][2])*u + c[i][1];
        ansy2 = t*ansy2 + (3.0f*c[i][4]*u + 2.0f*c[i][3])*u + c[i][2];
    }
    // For ansy1, note the polynomial in t is stored by columns in c.
    for (int j = 4; j >= 1; j--) {
        ansy1 = u*ansy1 + (3.0f*c[4][j]*t + 2.0f*c[3][j])*t + c[2][j];
    }

    // Adjust for the derivative scaling
    ansy1 /= d1;
    ansy2 /= d2;
}

int main() {
    srand((unsigned)time(nullptr));

    // We'll define a single cell from (x1l,y1l) = (0,0) to (x1u,y2u) = (1,1).
    float x1l = 0.0f, x1u = 1.0f;
    float x2l = 0.0f, x2u = 1.0f;

    // The corners, numbered counterclockwise from the lower left:
    //   1: (x1l, y1l) = (0,0)
    //   2: (x1u, y1l) = (1,0)
    //   3: (x1u, y2u) = (1,1)
    //   4: (x1l, y2u) = (0,1)
    //
    // We'll store function values f(x,y) and partial derivatives in arrays y[1..4], etc.

    float y[5], y1[5], y2[5], y12[5];

    // Corner 1: (0,0)
    y[1]   = fxy(0.0f, 0.0f);
    y1[1]  = cosf(0.0f)*cosf(0.0f);   // df/dx
    y2[1]  = -sinf(0.0f)*sinf(0.0f); // df/dy
    y12[1] = -cosf(0.0f)*sinf(0.0f); // d^2 f/(dx dy)

    // Corner 2: (1,0)
    y[2]   = fxy(1.0f, 0.0f);
    y1[2]  = cosf(1.0f)*cosf(0.0f);
    y2[2]  = -sinf(1.0f)*sinf(0.0f);
    y12[2] = -cosf(1.0f)*sinf(0.0f);

    // Corner 3: (1,1)
    y[3]   = fxy(1.0f, 1.0f);
    y1[3]  = cosf(1.0f)*cosf(1.0f);
    y2[3]  = -sinf(1.0f)*sinf(1.0f);
    y12[3] = -cosf(1.0f)*sinf(1.0f);

    // Corner 4: (0,1)
    y[4]   = fxy(0.0f, 1.0f);
    y1[4]  = cosf(0.0f)*cosf(1.0f);
    y2[4]  = -sinf(0.0f)*sinf(1.0f);
    y12[4] = -cosf(0.0f)*sinf(1.0f);

    // Now pick a random point (x, y) in [0,1] x [0,1] for interpolation
    float xTest = static_cast<float>(rand()) / RAND_MAX; // in [0,1]
    float yTest = static_cast<float>(rand()) / RAND_MAX; // in [0,1]

    // Call bcuint to interpolate at (xTest, yTest)
    bcuint(y, y1, y2, y12, x1l, x1u, x2l, x2u, xTest, yTest);

    // Compare with the true value
    float fTrue = fxy(xTest, yTest);

    // Print results
    cout << "Bicubic interpolation of f(x,y) = sin(x)*cos(y)\n\n";
    cout << "Grid corners: (0,0), (1,0), (1,1), (0,1)\n";
    cout << "Random point (xTest, yTest) = (" << xTest << ", " << yTest << ")\n\n";

    cout << "Interpolated value ansy   = " << ansy << "\n";
    cout << "True value f(xTest,yTest) = " << fTrue << "\n";
    cout << "Error                     = " << (ansy - fTrue) << "\n\n";

    cout << "Partial derivatives:\n";
    cout << "  df/dx (interpolated) = " << ansy1 << "\n";
    cout << "  df/dy (interpolated) = " << ansy2 << "\n";
    // For reference:
    float dfdx_true = cosf(xTest)*cosf(yTest);
    float dfdy_true = -sinf(xTest)*sinf(yTest);
    cout << "  df/dx (true)         = " << dfdx_true << "\n";
    cout << "  df/dy (true)         = " << dfdy_true << "\n";

    return 0;
}
