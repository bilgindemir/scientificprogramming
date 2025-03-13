//
// main.cpp
// sp_week1
//
// Created by Bilgin Demir on 2/27/25.
// Copyright (c) 2025 Bilgin Demir. All rights reserved.
//
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>

#include "Random.hpp"
#include "defs.hpp"
#include "Maths.hpp"

using namespace std;
using namespace rangen;
using namespace maths;

int main() {
    // Initialize random seed.
    init_random();

    int n;
    cout << "Enter the number of equations (n): ";
    cin >> n;
    while(cin.fail() || n <= 0) {
        cin.clear();
        cin.ignore(10000, '\n');
        cout << "Invalid input. Please enter a positive integer for n: ";
        cin >> n;
    }
    
    // Use SystemCount from defs.hpp.
    const int systemCount = SystemCount;

    // Ask the user to input the frequency (omega)
    float omega;
    cout << "Enter the frequency (omega): ";
    cin >> omega;
    
    // Generate a random time t in the range [0, 2*PI]
    float t = frand(0.0f, 2.0f * PI);
    cout << "Random time t generated: " << t << "\n" << endl;
    
    // Dynamically allocate matrices A (mass matrix) and b (force vector)
    // Using 1-indexing (allocate arrays of size n+1)
    float **A = new float*[n+1];
    float **b = new float*[n+1];
    for (int i = 0; i <= n; i++) {
        A[i] = new float[n+1];
        b[i] = new float[systemCount+1];
    }
    
    // Generate random mass matrix A:
    // A[i][j] = frand(0.0, 1.0) * MassMax
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            A[i][j] = frand(0.0f, 1.0f) * MassMax;
        }
    }
    
    // Generate the force vector b for each system:
    // b[i][k] = (frand(0.0,1.0)*2.0 + 1.0) * k * cos(omega * t)
    // Here k is the system index.
    for (int i = 1; i <= n; i++) {
        for (int k = 1; k <= systemCount; k++) {
            b[i][k] = (frand(0.0f, 1.0f)*2.0f + 1.0f) * k * cos(omega * t);
        }
    }
    
    // Display the generated mass matrix A.
    cout << "\nMass Matrix A:" << endl;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
    
    // Display the generated force vector b.
    cout << "\nForce Vector b:" << endl;
    for (int i = 1; i <= n; i++) {
        for (int k = 1; k <= systemCount; k++) {
            cout << b[i][k] << "\t";
        }
        cout << endl;
    }
    
    // Solve the system A * X = b to obtain the acceleration vector X.
    // The Gauss–Jordan elimination function overwrites b with the solution.
    gaussj(A, b, n, systemCount);
    
    // Display the computed acceleration vector X.
    cout << "\nAcceleration Vector X:" << endl;
    for (int i = 1; i <= n; i++) {
        for (int k = 1; k <= systemCount; k++) {
            cout << b[i][k] << "\t";
        }
        cout << endl;
    }
    
    // Free allocated memory.
    for (int i = 0; i <= n; i++) {
        delete[] A[i];
        delete[] b[i];
    }
    delete[] A;
    delete[] b;
    
    return 0;
}
