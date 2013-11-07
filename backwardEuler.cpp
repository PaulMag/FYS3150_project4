#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace arma;

vec backwardEuler() {

    int n = 100;
    double alpha = 0.123; // arbitrary value

    vec a(n);
    vec b(n);
    vec b_new(n);
    vec c(n);

    vec u     = zeros<vec>(n); // previous timestep
    vec u_new = zeros<vec>(n);
    vec v     = zeros<vec>(n); // next timestep
    u(0)     = 1;
    u_new(0) = 1;
    v(0)     = 1;

    // Fill in arrays:
    for (int i=0; i<n; i++) {
        a(i) = -alpha;
        b(i) = 1 + 2 * alpha;
        c(i) = -alpha;
    }
    a(0)   = 0; // parts that "stick outside" matrix
    c(n-1) = 0;

    // The algorithm:
    b_new(0) = b(0);
    u_new(0) = u(0);

    for (int i=1; i<n; i++) {
        double factor = a(i) / b_new(i-1); // avoids doing this twice

        b_new(i) = b(i) - c(i-1) * factor;

        u_new(i) = u(i) - u_new(i-1) * factor;
    }

    v(n-1) = u(n-1) / b_new(n-1);
    for (int i=n-2; i>0; i--) {
        v(i) = (u_new(i) - v(i+1) * c(1)) / b_new(i); //unsure which c,no matter
    }

    return v;
}
