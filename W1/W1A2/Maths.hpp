#ifndef MATHS_HPP
#define MATHS_HPP

#include <iostream>
#include <cmath>

namespace maths {

// Macro for swapping two values
#define SWAP(a,b) { float temp = (a); (a) = (b); (b) = temp; }

// Gauss–Jordan elimination function.
// Solves A * X = b for X. Both matrices are assumed to be 1-indexed.
// On output, A is replaced by its inverse and b by the solution.
void gaussj(float **a, float **b, int n, int m) {
    int i, icol, irow, j, k, l, ll;
    float big, dum, pivinv;

    // Allocate temporary arrays (1-indexed)
    int *indxc = new int[n+1];
    int *indxr = new int[n+1];
    int *ipiv  = new int[n+1];

    for (j = 1; j <= n; j++)
        ipiv[j] = 0;

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

        // Swap rows if needed.
        if (irow != icol) {
            for (l = 1; l <= n; l++)
                SWAP(a[irow][l], a[icol][l]);
            for (l = 1; l <= m; l++)
                SWAP(b[irow][l], b[icol][l]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) {
            std::cerr << "gaussj: Singular Matrix" << std::endl;
            delete[] indxc; delete[] indxr; delete[] ipiv;
            return;
        }
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++)
            a[icol][l] *= pivinv;
        for (l = 1; l <= m; l++)
            b[icol][l] *= pivinv;
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

} // namespace maths

#endif // MATHS_HPP
