//
// main.cpp
// assignment1_week2
//
// Created by Hiqmet Kamberaj on 4/13/23.
// Copyright (c) 2023 Hiqmet Kamberaj. All rights reserved.
// 

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>

// Define the function we want to interpolate/extrapolate

float fun(float x) {
    return 3.0f * x * x
         + (1.0f / (float)std::pow(M_PI, 4)) * std::log((M_PI - x)*(M_PI - x))
         + 1.0f;
}

// A simple selection sort to keep our x-array in ascending order
// (Numerical Recipes typically uses 1-based indexing, so we do too.)

void selectionSort(int start, int end, float arr[]) {
    for (int i = start; i < end; i++) {
        int minIndex = i;
        for (int j = i + 1; j <= end; j++) {
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }
        // Swap
        float temp = arr[i];
        arr[i] = arr[minIndex];
        arr[minIndex] = temp;
    }
}

// The Numerical Recipes polint function (polynomial interpolation)

float y, dy;  // Global variables, as per Numerical Recipes convention

void polint(float xa[], float ya[], int n, float x) {
    int i, m, ns = 1;
    float den, dif, dift, ho, hp, w;

    // Allocate local workspace for c and d
    float *c = new float[n+1];
    float *d = new float[n+1];

    // Find the closest point in xa[] to x, initialize c[] and d[]
    dif = std::fabs(x - xa[1]);
    for (i = 1; i <= n; i++) {
        dift = std::fabs(x - xa[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    // Initial approximation
    y = ya[ns--];

    // Perform the nested interpolation
    for (m = 1; m < n; m++) {
        for (i = 1; i <= n - m; i++) {
            ho = xa[i]     - x;
            hp = xa[i + m] - x;
            w  = c[i + 1]  - d[i];
            den = ho - hp;
            if (den == 0.0f) {
                std::cerr << "Error in polint: den=0 (identical x-values?)\n";
                delete[] c;
                delete[] d;
                return;
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        dy = (2 * ns < (n - m)) ? c[ns + 1] : d[ns--];
        y += dy;
    }

    delete[] c;
    delete[] d;
}

// Main program:
//    - Generate n random points (x_i, y_i) in ascending order
//    - Let the user pick interpolation vs. extrapolation
//    - Compute polint(x0) and compare to the "true" f(x0)
int main() {
    // Initialize random seed
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    // Number of data points
    int n = 10;

    // 1-based arrays to store data
    float xa[n+1], ya[n+1];

    // Generate random x-values in [0, 0.5 / π], for example
    for (int i = 1; i <= n; i++) {
        xa[i] = (float)std::rand() / RAND_MAX;   // in [0,1]
        xa[i] *= 0.5f * (1.0f / M_PI);           // scale by 0.5 / π
    }

    // Sort the x-values
    selectionSort(1, n, xa);

    // Compute f(x) at these points
    for (int i = 1; i <= n; i++) {
        ya[i] = fun(xa[i]);
    }

    // Print known data points
    std::cout << "Known data points (x_i, y_i):\n";
    for (int i = 1; i <= n; i++) {
        std::cout << "  x[" << i << "] = " << xa[i]
                  << ", y[" << i << "] = " << ya[i] << "\n";
    }

    // Ask user whether to interpolate or extrapolate
    int choice;
    std::cout << "\nChoose:\n"
              << "  1) Interpolation (x0 inside [x1, xN])\n"
              << "  2) Extrapolation (x0 outside [x1, xN])\n"
              << "Enter choice: ";
    std::cin >> choice;
    while (std::cin.fail() || (choice != 1 && choice != 2)) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "Invalid input. Please enter 1 or 2: ";
        std::cin >> choice;
    }

    // Determine the min and max x-values
    float xMin = xa[1];
    float xMax = xa[n];

    // Decide on x0 for interpolation/extrapolation
    float x0;
    if (choice == 1) {
        // Interpolation: pick x0 in [xMin, xMax]
        x0 = xMin + ((float)std::rand() / RAND_MAX) * (xMax - xMin);
        std::cout << "\nInterpolation selected. x0 in [" 
                  << xMin << ", " << xMax << "]\n";
    } else {
        // Extrapolation: pick x0 below xMin or above xMax
        float belowOrAbove = (float)std::rand() / RAND_MAX;
        if (belowOrAbove < 0.5f) {
            // pick below xMin
            x0 = xMin - 0.1f * (xMax - xMin);
            std::cout << "\nExtrapolation selected. x0 chosen below " 
                      << xMin << ".\n";
        } else {
            // pick above xMax
            x0 = xMax + 0.1f * (xMax - xMin);
            std::cout << "\nExtrapolation selected. x0 chosen above " 
                      << xMax << ".\n";
        }
    }

    std::cout << "x0 = " << x0 << "\n";

    // Interpolate (or extrapolate) using polint
    polint(xa, ya, n, x0);

    // Compare to the "true" function value
    float trueVal = fun(x0);

    // Print results
    std::cout << "\nResults:\n"
              << "  x0 = " << x0 << "\n"
              << "  True f(x0)  = " << trueVal << "\n"
              << "  polint(x0)  = " << y << " +/- " << dy << "\n\n";

    return 0;
}
