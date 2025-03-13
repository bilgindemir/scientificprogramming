///////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Author: Bilgin Demir
 *
 * Created on 27 February 2025
 */
/////////////////////////////////////////////////////////////////////////////////////////////////// 
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>

using namespace std;

// Macro for swapping two values.
#define SWAP(a,b) { float temp = (a); (a) = (b); (b) = temp; }

// Helper function to print a matrix (for ASCII visualization)
void printMatrix(float **matrix, int rows, int cols) {
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
}

// Gauss-Jordan elimination function.
// Solves A x = b by replacing A with its inverse and b with the solution vector.
// 'a' is the coefficient matrix (1-indexed) and 'b' is the right-hand side matrix.
void gaussj(float **a, float **b, int n, int m) {
    int i, icol, irow, j, k, l, ll;
    float big, dum, pivinv;

    // Allocate temporary arrays for bookkeeping (1-indexed)
    int *indxc = new int[n+1];
    int *indxr = new int[n+1];
    int *ipiv  = new int[n+1];

    // Initialize pivot markers.
    for (j = 1; j <= n; j++) 
        ipiv[j] = 0;

    // Main loop over the columns.
    for (i = 1; i <= n; i++) {
        big = 0.0;
        // Search for the pivot element.
        for (j = 1; j <= n; j++) {
            if (ipiv[j] != 1) {
                for (k = 1; k <= n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ipiv[icol]++;

        // Swap rows if necessary.
        if (irow != icol) {
            for (l = 1; l <= n; l++)
                SWAP(a[irow][l], a[icol][l]);
            for (l = 1; l <= m; l++)
                SWAP(b[irow][l], b[icol][l]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) {
            cout << "gaussj: Singular Matrix" << endl;
            delete[] indxc; delete[] indxr; delete[] ipiv;
            return;
        }
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++)
            a[icol][l] *= pivinv;
        for (l = 1; l <= m; l++)
            b[icol][l] *= pivinv;

        // Eliminate the pivot column from the other rows.
        for (ll = 1; ll <= n; ll++) {
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 1; l <= n; l++)
                    a[ll][l] -= a[icol][l] * dum;
                for (l = 1; l <= m; l++)
                    b[ll][l] -= b[icol][l] * dum;
            }
        }
    }
    
    // Undo column interchanges.
    for (l = n; l >= 1; l--) {
        if (indxr[l] != indxc[l]) {
            for (k = 1; k <= n; k++)
                SWAP(a[k][indxr[l]], a[k][indxc[l]]);
        }
    }
    delete[] indxc;
    delete[] indxr;
    delete[] ipiv;
}

int main() {
    int n, m = 1;  // m = 1 because we have one right-hand side vector

    // Prompt for the number of equations (and unknowns).
    cout << "Enter the number of equations (n), which is also the number of unknowns: ";
    cin >> n;
    
    // Validate input to ensure n is a positive integer.
    while (cin.fail() || n <= 0) {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Invalid input. Please enter a positive integer for the number of equations: ";
        cin >> n;
    }
    
    // Dynamically allocate matrices 'a' (n x n) and 'b' (n x m) using 1-indexing.
    float **a = new float*[n+1];
    float **b = new float*[n+1];
    for (int i = 0; i <= n; i++) {
        a[i] = new float[n+1];
        b[i] = new float[m+1];
    }
    
    // Read matrix A from user input.
    cout << "\nEnter the coefficients of matrix A (row-wise):\n";
    for (int i = 1; i <= n; i++) {
        cout << "Row " << i << " (enter " << n << " coefficients separated by spaces): ";
        for (int j = 1; j <= n; j++) {
            cin >> a[i][j];
        }
    }
    
    // Read the right-hand side vector b.
    cout << "\nEnter the right-hand side vector b:\n";
    for (int i = 1; i <= n; i++) {
        cout << "b[" << i << "]: ";
        cin >> b[i][1];
    }
    
    // Display the input matrices (for visualization/debugging).
    cout << "\nMatrix A:" << endl;
    printMatrix(a, n, n);
    cout << "\nVector b:" << endl;
    for (int i = 1; i <= n; i++) {
        cout << b[i][1] << "\t";
    }
    cout << "\n\n";
    
    // Solve the linear system A x = b using Gauss-Jordan elimination.
    gaussj(a, b, n, m);
    
    // Output the solution vector.
    cout << "Solution vector x:" << endl;
    for (int i = 1; i <= n; i++) {
        cout << "x[" << i << "] = " << b[i][1] << endl;
    }
    
    // Free dynamically allocated memory.
    for (int i = 0; i <= n; i++) {
        delete[] a[i];
        delete[] b[i];
    }
    delete[] a;
    delete[] b;
    
    return 0;
}