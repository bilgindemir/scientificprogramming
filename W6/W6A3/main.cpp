#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

// Returns ln[Gamma(xx)] for xx > 0 via Lanczos approximation
double gammln(double xx) {
    static double cof[6] = {
        76.18009172947146,  -86.50532032941677,
        24.01409824083091,  -1.231739572450155,
         0.1208650973866179e-2,  -0.5395239384953e-5
    };
    double x = xx, y = xx;
    double tmp = x + 5.5;
    tmp -= (x + 0.5)*log(tmp);
    double ser = 1.000000000190015;
    for (int j = 0; j < 6; j++)
        ser += cof[j]/(++y);
    return -tmp + log(2.5066282746310005 * ser / x);
}

// Returns n! exactly for n<=32, else via exp(gammln(n+1))
double factorial(int n) {
    static int   ntop = 4;
    static double a[33] = {
        1.0,   1.0,    2.0,     6.0,
       24.0, 120.0, 720.0, 5040.0,
    // remaining entries will be filled on demand
    };
    if (n < 0) {
        cerr << "Negative factorial\n";
        return 1.0;
    }
    if (n > 32) {
        return exp(gammln(n+1.0));
    }
    // fill table up to n
    while (ntop < n) {
        int j = ntop++;
        a[ntop] = a[j]*ntop;
    }
    return a[n];
}

// Returns ln(n!). Caches values up to n=100.
double factorialln(int n) {
    static double a[101] = {0.0};
    if (n < 0) {
        cerr << "Negative factorialln\n";
        return 0.0;
    }
    if (n <= 1) return 0.0;
    if (n <= 100) {
        if (a[n] != 0.0) return a[n];
        return a[n] = gammln(n+1.0);
    }
    return gammln(n+1.0);
}

// Returns binomial coefficient C(n,k) as a double
double bico(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    double lnC = factorialln(n) - factorialln(k) - factorialln(n - k);
    return floor(0.5 + exp(lnC));
}

// Returns Beta(z,w) = Gamma(z)*Gamma(w)/Gamma(z+w)
double beta(double z, double w) {
    return exp( gammln(z) + gammln(w) - gammln(z+w) );
}

int main(){
    cout << fixed << setprecision(9);

    // 1) ln(Gamma)
    cout << "ln Gamma(5.5) = " << gammln(5.5)
         << "  (exact ~ " << log( gamma(5.5) ) << ")\n";

    // 2) factorial
    cout << "5! = " << factorial(5) << "\n";
    cout << "33! (via gammln) = " << factorial(33) << "\n";

    // 3) factorialln
    cout << "ln(10!) = " << factorialln(10)
         << "  (exact ~ " << log(3628800.0) << ")\n";

    // 4) binomial coefficients
    cout << "C(10,3) = " << bico(10,3) << "  (should be 120)\n";
    cout << "C(50,25) = " << bico(50,25) << "\n";

    // 5) Beta function
    cout << "Beta(2.5,3.5) = " << beta(2.5,3.5)
         << "  (exact ~ " << tgamma(2.5)*tgamma(3.5)/tgamma(6.0) << ")\n";

    return 0;
}
