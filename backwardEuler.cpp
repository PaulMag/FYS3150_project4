#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace arma;

vec backwardEuler() {

    int n = 100;
    double alpha = 0.123; // arbitrary value

    vec a(n-1);
    vec b(n-1);     // previous timestep
    vec b_new(n-1);
    vec c(n-1);     // next timestep

    vec u     = zeros<vec>(n+1);
    vec u_new = zeros<vec>(n+1);
    vec v     = zeros<vec>(n+1);
    u(0)     = 1;
    u_new(0) = 1;
    v(0)     = 1;

    // Fill in arrays:
    for (int i=0; i<n-1; i++) {
        a(i) = -alpha;
        b(i) = 1 + 2 * alpha;
        c(i) = -alpha;
    }
    a(0)   = 0; // parts that "stick outside" matrix
    c(n-2) = 0;

    // The algorithm:
    for (int i=1; i<n-1; i++) {
        double factor = a(i) / b_new(i-1); // avoids doing this twice

        b_new(i) = b(i) - c(i-1) * factor;

        u_new(i +1) = u(i +1) - u_new(i-1 +1) * factor;
    }

    //v(n-1) = u(n-1) / b_new(n-1);
    for (int i=n-2; i>=0; i--) {
        v(i +1) = (u_new(i +1) - v(i+1 +1) * c(i)) / b_new(i);
    }

    return v;
}
