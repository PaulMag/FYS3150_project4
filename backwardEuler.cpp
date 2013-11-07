#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace arma;

vec backwardEuler(int n, double time, double dt) {

    double alpha = dt * (n+1) * (n+1); // alpha = dt / dx^2

    vec a(n);
    vec b(n);
    vec b_new(n);
    vec c(n);

    // Fill in arrays:
    for (int i=0; i<n; i++) {
        a(i) = -alpha;
        b(i) = 1 + 2 * alpha;
        c(i) = -alpha;
    }
    a(0)   = 0; // parts that "stick outside" matrix
    c(n-1) = 0; // no point in this a.t.m.

    // Initial value:
    vec u = zeros<vec>(n); // previous timestep
    u(0)  = 1;
    vec u_new(n);
    vec v(n); // next timestep
    v(0)  = u(0);

    cout << u << endl; // for testing

    // Time loop:
    for (double j=0; j*dt < time; j++) {

        // Algorithm:
        b_new(0) = b(0);
        u_new(0) = u(0);

        for (int i=1; i<n; i++) {
            double factor = a(i) / b_new(i-1); // avoids doing this twice

            b_new(i) = b(i) - c(i-1) * factor;

            u_new(i) = u(i) - u_new(i-1) * factor;
        }

        v(n-1) = u(n-1) / b_new(n-1);
        for (int i=n-2; i>0; i--) {
            v(i) = (u_new(i) - v(i+1) * c(i)) / b_new(i); //unsure which c,no matter
        }

        for (int i=0; i<n; i++) {
            u(i) = v(i);
        }
        //cout << u << endl; // for testing
    }

    return v;
}
