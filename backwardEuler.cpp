#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace arma;

vec backwardEuler(int n, double time, double dt) {
    /* @param n Resolution in position: n = 1/dt
     * @param time How long to integrate.
     * @param dt Timestep.
     *
     * The _alt extension on a vec means that it is altered through row
     * reduction. u_new( is how u looks in the next timestep, what we want to
     * find.
     */

    double alpha = dt * (n+1) * (n+1); // alpha = dt / dx^2

    vec a(n);
    vec b(n);
    vec b_alt(n);
    vec c(n);

    // Fill in arrays:
    for (int i=0; i<n; i++) {
        a(i) = -alpha;
        b(i) = 1 + 2 * alpha;
        c(i) = -alpha;
    }

    // Initial value:
    vec u = zeros<vec>(n); // previous timestep
    double u0  = 1.0;      // this is the element before u(0)
    vec u_alt(n);
    vec u_new(n); // next timestep


    // Time loop:
    for (double j=0; j*dt < time; j++) {

        // Algorithm:
        b_alt(0) = b(0);
        u_alt(0) = u(0) - u0 * a(0);

        for (int i=1; i<n; i++) {
            double factor = a(i) / b_alt(i-1); // avoids doing this twice

            b_alt(i) = b(i) - factor * c(i-1);
            u_alt(i) = u(i) - factor * u_alt(i-1);
        }

        u_new(n-1) = u(n-1) / b_alt(n-1);
        for (int i=n-2; i>=0; i--) {
            u_new(i) = (u_alt(i) - u_new(i+1) * c(i)) / b_alt(i);
        }

        for (int i=0; i<n; i++) {
            u(i) = u_new(i);
        }
    }

    return u_new;
}
