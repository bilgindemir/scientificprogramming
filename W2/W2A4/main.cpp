#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <limits>

using namespace std;

// -----------------------------------------------------------------------------
// 1) The function to be interpolated
//    f(x) = 3x^2 + (1 / π^4) * log((π - x)^2) + 1
// -----------------------------------------------------------------------------
float fun(float x) {
    return 3.0f * x * x
         + (1.0f / (float)pow(M_PI, 4)) * log((M_PI - x)*(M_PI - x))
         + 1.0f;
}

// -----------------------------------------------------------------------------
// 2) Evaluate a polynomial with coefficients cof[0..n] at point x
//    P(x) = cof[0] + cof[1]*x + cof[2]*x^2 + ... + cof[n]*x^n
// -----------------------------------------------------------------------------
float evaluatePoly(const vector<float>& cof, float x) {
    float val = 0.0f;
    float xn  = 1.0f;
    for (int i = 0; i < (int)cof.size(); i++) {
        val += cof[i] * xn;
        xn  *= x;
    }
    return val;
}

// -----------------------------------------------------------------------------
// 3) polcoe1: Direct Lagrange-type polynomial coefficients
//    Arrays x[0..n], y[0..n] => cof[0..n]
//    Based on Numerical Recipes code snippet (polcoe1).
// -----------------------------------------------------------------------------
void polcoe1(const vector<float>& xv, const vector<float>& yv, int n, vector<float>& cof) {
    // Make sure cof is sized correctly
    cof.assign(n+1, 0.0f);

    // Temporary "master polynomial" array
    vector<float> s(n+1, 0.0f);

    // Initialize
    for (int i = 0; i <= n; i++) {
        s[i] = 0.0f;
        cof[i] = 0.0f;
    }

    // Build the "master polynomial" s[n..0]
    s[n] = -xv[0];
    for (int i = 1; i <= n; i++) {
        for (int j = n - i; j <= n - 1; j++) {
            s[j] -= xv[i] * s[j + 1];
        }
        s[n] -= xv[i];
    }

    // Construct the polynomial coefficients
    for (int j = 0; j <= n; j++) {
        float phi = n + 1;
        for (int k = n; k >= 1; k--) {
            phi = k * s[k] + xv[j] * phi;
        }
        float ff = yv[j] / phi;

        float b = 1.0f;
        for (int k = n; k >= 0; k--) {
            cof[k] += b * ff;
            b = s[k] + xv[j] * b;
        }
    }
}

// -----------------------------------------------------------------------------
// 4) polcoe2: Iterative approach that uses polint to find each coefficient
//    by extrapolating to x=0, then modifying the data, repeating.
//    This snippet references a polint routine with 1-based indexing.
// -----------------------------------------------------------------------------

// polint adapted for 0-based indexing here for convenience.
void polint(const float* xa, const float* ya, int n, float x, float& y, float& dy) {
    vector<float> c(n), d(n);

    // Initialize
    int ns = 0;
    float dif = fabs(x - xa[0]);
    for (int i = 0; i < n; i++) {
        float dift = fabs(x - xa[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    y = ya[ns];
    ns--;

    for (int m = 1; m < n; m++) {
        for (int i = 0; i < n - m; i++) {
            float ho = xa[i]     - x;
            float hp = xa[i + m] - x;
            float w  = c[i + 1]  - d[i];
            float den = ho - hp;
            if (den == 0.0f) {
                cerr << "Error in polint: den=0\n";
                return;
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        if ((2 * ns + 2) < (n - m)) {
            dy = c[ns + 1];
        } else {
            dy = d[ns];
            ns--;
        }
        y += dy;
    }
}

void polcoe2(const vector<float>& xa, const vector<float>& ya, int n, vector<float>& cof) {
    // We'll copy the data since polcoe2 modifies it
    vector<float> xx = xa;
    vector<float> yy = ya;

    cof.assign(n+1, 0.0f);

    for (int j = 0; j <= n; j++) {
        // polint expects 0-based indexing, so pass addresses
        float c0, dy;
        polint(&xx[0], &yy[0], n + 1 - j, 0.0f, c0, dy);
        cof[j] = c0;

        // Subtract c0 from y-values, divide by x[i], remove one point
        // pick the point with the smallest |x[i]| to remove
        float xmin = 1.0e38f;
        int k = -1;
        for (int i = 0; i <= n - j; i++) {
            if (fabs(xx[i]) < xmin) {
                xmin = fabs(xx[i]);
                k = i;
            }
            if (xx[i] != 0.0f) {
                yy[i] = (yy[i] - c0) / xx[i];
            }
        }
        // "Remove" that point
        for (int i = k + 1; i <= n - j; i++) {
            xx[i - 1] = xx[i];
            yy[i - 1] = yy[i];
        }
    }
}

// -----------------------------------------------------------------------------
// 5) Main test program
//    - Generates n random data points in [0,1], computes f(x).
//    - Gets polynomial coefficients by polcoe1 and polcoe2.
//    - Evaluates both polynomials at random test points, compares with f(x).
//    - Prints max error and RMS error for each method.
// -----------------------------------------------------------------------------
int main() {
    srand(static_cast<unsigned>(time(nullptr)));

    // Choose how many data points
    int n = 20;  // You can try bigger n, e.g., 15 or 20, to see more numerical issues

    // Generate x in [0,1] (you can try other intervals)
    vector<float> x(n+1), y(n+1);
    for (int i = 0; i <= n; i++) {
        x[i] = static_cast<float>(rand()) / RAND_MAX; // in [0,1]
        y[i] = fun(x[i]);
    }

    // polcoe1 vs. polcoe2
    vector<float> cof1, cof2;
    polcoe1(x, y, n, cof1);
    polcoe2(x, y, n, cof2);

    // Evaluate at some test points
    int nTest = 50; // number of test points
    double maxErr1 = 0.0, maxErr2 = 0.0;
    double rmsErr1 = 0.0, rmsErr2 = 0.0;

    for (int i = 0; i < nTest; i++) {
        // random test point in [0,1]
        float xt = static_cast<float>(rand()) / RAND_MAX;
        float ftrue = fun(xt);

        float f1 = evaluatePoly(cof1, xt);
        float f2 = evaluatePoly(cof2, xt);

        double err1 = fabs(f1 - ftrue);
        double err2 = fabs(f2 - ftrue);

        if (err1 > maxErr1) maxErr1 = err1;
        if (err2 > maxErr2) maxErr2 = err2;

        rmsErr1 += err1 * err1;
        rmsErr2 += err2 * err2;
    }

    rmsErr1 = sqrt(rmsErr1 / nTest);
    rmsErr2 = sqrt(rmsErr2 / nTest);

    cout << "\n=== Polynomial Coefficients Stability Test ===\n";
    cout << "Data points: n = " << n+1 << "\n";
    cout << "Test points: " << nTest << "\n\n";

    cout << "Method 1 (polcoe1):\n";
    cout << "  Max error = " << maxErr1 << "\n";
    cout << "  RMS error = " << rmsErr1 << "\n\n";

    cout << "Method 2 (polcoe2):\n";
    cout << "  Max error = " << maxErr2 << "\n";
    cout << "  RMS error = " << rmsErr2 << "\n\n";

    if (maxErr2 > maxErr1 * 10.0) {
        cout << "=> polcoe2 shows higher error. This illustrates lower stability.\n";
    } else {
        cout << "=> In this run, polcoe2 error is not drastically worse. (Try larger n or different x distribution.)\n";
    }

    return 0;
}
