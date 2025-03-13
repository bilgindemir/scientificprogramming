//
// main.cpp
// assignment2_week2
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

//--------------------------------------------------------------------------
// Global constant to prevent a zero denominator in ratint.
const float TINY = 1.0e-25f;

//--------------------------------------------------------------------------
// Define the function to be tabulated.
// For example, we use the function from Assignment 1:
//    f(x) = 3x^2 + (1/π^4)*log((π - x)^2) + 1
// (You can change this to any function for which you have tabulated data.)
float fun(float x) {
    return 3.0f * x * x +
           (1.0f / (float)pow(M_PI, 4)) * log((M_PI - x) * (M_PI - x)) +
           1.0f;
}

//--------------------------------------------------------------------------
// Simple selection sort to sort the array of x values (1-indexed).
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

//--------------------------------------------------------------------------
// Global variables to store the result and error estimate
// (following the Numerical Recipes convention).
float y, dy;

//--------------------------------------------------------------------------
// Rational interpolation routine (ratint).
// Given arrays xa[1..n] and ya[1..n], and a value x, this routine
// returns an interpolated (or extrapolated) value y and an error estimate dy.
// The diagonal rational function constructed passes through the n points.
void ratint(float xa[], float ya[], int n, float x) {
    int m, i, ns = 1;
    float w, t, hh, h, dd;
    // Create local arrays for c and d (1-indexed)
    float *c = new float[n+1];
    float *d = new float[n+1];
    
    hh = fabs(x - xa[1]);
    for (i = 1; i <= n; i++) {
        h = fabs(x - xa[i]);
        if (h == 0.0f) {
            y = ya[i];
            dy = 0.0f;
            delete[] c;
            delete[] d;
            return;
        } else if (h < hh) {
            ns = i;
            hh = h;
        }
        c[i] = ya[i];
        d[i] = ya[i] + TINY;  // Prevent a zero-over-zero condition.
    }
    y = ya[ns--];
    for (m = 1; m < n; m++) {
        for (i = 1; i <= n - m; i++) {
            w = c[i+1] - d[i];
            // Here h should not be zero because it was tested above.
            h = xa[i+m] - x;
            t = (xa[i] - x) * d[i] / h;
            dd = t - c[i+1];
            if (dd == 0.0f) {
                cout << "Error in routine ratint: Pole encountered at x = " << x << endl;
                // We continue; the error condition signals a pole.
            }
            dd = w / dd;
            d[i] = c[i+1] * dd;
            c[i] = t * dd;
        }
        dy = (2 * ns < (n - m)) ? c[ns+1] : d[ns--];
        y += dy;
    }
    delete[] c;
    delete[] d;
}

//--------------------------------------------------------------------------
// Main program
// - Generates n data points (xi, yi) with xi sorted in ascending order.
// - Prompts the user to choose between interpolation (x0 inside [x1, xN])
//   or extrapolation (x0 outside [x1, xN]).
// - Uses ratint to estimate f(x0) and compares it to the true value.
int main() {
    // Initialize random seed.
    srand(static_cast<unsigned>(time(nullptr)));
    
    const int n = 10;  // Number of tabulated data points.
    float xa[n+1], ya[n+1];
    
    // Generate n random x-values in an interval.
    // For this example, we choose x in [0, 0.5 / π]
    for (int i = 1; i <= n; i++) {
        xa[i] = static_cast<float>(rand()) / RAND_MAX; // in [0,1]
        xa[i] *= 0.5f * (1.0f / M_PI);                  // scale to [0, 0.5/π]
    }
    
    // Sort the x-values (using 1-indexing).
    selectionSort(1, n, xa);
    
    // Compute corresponding y-values using the function fun.
    for (int i = 1; i <= n; i++) {
        ya[i] = fun(xa[i]);
    }
    
    // Print out the known data points.
    cout << "Tabulated data points (xi, f(xi)):" << endl;
    for (int i = 1; i <= n; i++) {
        cout << "  x[" << i << "] = " << xa[i]
             << "   f(x[" << i << "]) = " << ya[i] << endl;
    }
    
    // Ask the user whether to perform interpolation or extrapolation.
    int choice;
    cout << "\nChoose an option:" << endl;
    cout << "  1) Interpolation (choose x0 inside the range [x1, xN])" << endl;
    cout << "  2) Extrapolation (choose x0 outside the range [x1, xN])" << endl;
    cout << "Enter choice (1 or 2): ";
    cin >> choice;
    while (cin.fail() || (choice != 1 && choice != 2)) {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Invalid input. Please enter 1 or 2: ";
        cin >> choice;
    }
    
    // Determine the minimum and maximum x values from the data.
    float xMin = xa[1];
    float xMax = xa[n];
    float x0;
    if (choice == 1) {
        // Interpolation: choose x0 inside [xMin, xMax]
        x0 = xMin + (static_cast<float>(rand()) / RAND_MAX) * (xMax - xMin);
        cout << "\nInterpolation selected. x0 is chosen within [" 
             << xMin << ", " << xMax << "]." << endl;
    } else {
        // Extrapolation: choose x0 outside the interval.
        // Here, we pick x0 either a bit below xMin or above xMax.
        float r = static_cast<float>(rand()) / RAND_MAX;
        if (r < 0.5f) {
            x0 = xMin - 0.1f * (xMax - xMin);
            cout << "\nExtrapolation selected. x0 is chosen below " << xMin << "." << endl;
        } else {
            x0 = xMax + 0.1f * (xMax - xMin);
            cout << "\nExtrapolation selected. x0 is chosen above " << xMax << "." << endl;
        }
    }
    cout << "x0 = " << x0 << endl;
    
    // Perform rational interpolation (or extrapolation) at x0.
    ratint(xa, ya, n, x0);
    
    // Compute the "true" function value at x0.
    float trueVal = fun(x0);
    
    // Output the results.
    cout << "\nResults:" << endl;
    cout << "  x0 = " << x0 << endl;
    cout << "  True f(x0)  = " << trueVal << endl;
    cout << "  ratint(x0)  = " << y << " +/- " << dy << endl;
    
    return 0;
}
