#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace arma;


rowvec analytic(int n, double time) {

    int N = 10000;

    rowvec u(n+1);
    rowvec x(n+1);

    double dx = 1. / n;
    for (int i=0; i<n+1; i++) { // linspace
        x(i) = i * dx;
    }

    double pi  = 3.14159265359;
    double pi2 = pi * pi;

    for (int i=0; i<n; i++) {

        double s = 0.0;

        for (int k=1; k<=N; k++) {
            s += sin(k*pi*x(i)) / k * exp(- k*k * pi2 * time);
        }
        u(i) =  1 - x(i) - 2./pi * s;
    }

    return u;
}


rowvec forwardEuler(int n, double time, double dt) {

    double alpha = dt * n * n; // alpha = dt / dx^2
    if (alpha > 0.5) {
        cout << "WARNING: alpha = " << alpha << " < 1/2" << endl
             << "WARNING: Algorithm may be unstable!"    << endl;
    }

    // Initial value:
    rowvec u = zeros<rowvec>(n+1); // previous timestep
    u(0)  = 1.0;

    rowvec u_new(n+1); // next timestep
    u_new(0) = u(0);
    u_new(n) = u(n);

    // Time loop:
    for (double j=0; j*dt < time; j++) {

        // Algorithm:
        for (int i=1; i<n; i++) {
            u_new(i) = u(i) + ( u(i+1) - 2 * u(i) + u(i-1) ) * alpha;
        }

        // Update:
        u = u_new;
    }

    return u;
}


rowvec backwardEuler(int n, double time, double dt) {
    /* @param n Resolution in position: n = 1/dt
     * @param time How long to integrate.
     * @param dt Timestep.
     *
     * The _alt extension on a rowvec means that it is altered through row
     * reduction. u_new( is how u looks in the next timestep, what we want to
     * find.
     */

    double alpha = dt * n * n; // alpha = dt / dx^2

    rowvec a(n);
    rowvec b(n);
    rowvec b_alt(n);
    rowvec c(n);

    // Fill in arrays:
    for (int i=0; i<n; i++) {
        a(i) = -alpha;
        b(i) = 1 + 2 * alpha;
        c(i) = -alpha;
    }

    // Initial value:
    rowvec u     = zeros<rowvec>(n); // previous timestep
    double u0 = 1.0; // this is "u(-1)" the element before u(0)
    rowvec u_alt(n);
    rowvec u_new(n); // next timestep


    // Time loop:
    for (double j=0; j*dt < time; j++) {

        // Algorithm:
        b_alt(0) = b(0);
        u_alt(0) = u(0) - u0 * a(0);

        for (int i=1; i<n; i++) {
            /* Forwards substitution. */
            double factor = a(i) / b_alt(i-1); // avoids doing this twice

            b_alt(i) = b(i) - factor * c(i-1);
            u_alt(i) = u(i) - factor * u_alt(i-1);
        }

        u_new(n-1) = u(n-1) / b_alt(n-1);
        for (int i=n-2; i>=0; i--) {
            /* Backwards substitution. */
            u_new(i) = (u_alt(i) - u_new(i+1) * c(i)) / b_alt(i);
        }

        // Update:
        u = u_new;
    }

    // Add u0 at the beginning of u:
    rowvec u_fin(n+1);
    for (int i=0; i<n; i++) {
        u_fin(i+1) = u(i);
    }
    u_fin(0) = u0;

    return u_fin;
}


rowvec crankNicolson(int n, double time, double dt) {
    /* Takes half a step with forwardEuler and another half with backwardEuler.
     * WARNING: Length of the returned vec is only n, because the u0 at the
     * beginning is missing. This was the only way to fix the weirdest bug in
     * the universe.
     */

    double alpha = 0.5 * dt * n * n; // this is alpha/2

    rowvec a(n);
    rowvec b(n);
    rowvec b_alt(n);
    rowvec c(n);

    // Fill in arrays:
    for (int i=0; i<n; i++) {
        a(i) = -alpha;
        b(i) = 1 + 2 * alpha;
        c(i) = -alpha;
    }

    // Initial value:
    rowvec u     = zeros<rowvec>(n); // previous timestep
    double u0 = 1.0; // this is "u(-1)" the element before u(0)
    rowvec u_alt(n);
    rowvec u_mid(n); // half timestep
    rowvec u_new(n); // next timestep


    // Time loop:
    for (double j=0; j*dt < time; j++) {

        // Forward euler algorithm:
        u_mid(0) = u(0) + ( u(1) - 2 * u(0) + u0 ) * alpha; // first element
        for (int i=1; i<n-1; i++) {
            u_mid(i) = u(i) + ( u(i+1) - 2 * u(i) + u(i-1) ) * alpha;
        }

        // Backward euler algorithm:
        b_alt(0) = b(0);
        u_alt(0) = u_mid(0) - u0 * a(0);

        for (int i=1; i<n; i++) {
            /* Forwards substitution. */
            double factor = a(i) / b_alt(i-1); // avoids doing this twice

            b_alt(i) = b(i)     - factor * c(i-1);
            u_alt(i) = u_mid(i) - factor * u_alt(i-1);
        }

        u_new(n-1) = u_mid(n-1) / b_alt(n-1);
        for (int i=n-2; i>=0; i--) {
            /* Backwards substitution. */
            u_new(i) = (u_alt(i) - u_new(i+1) * c(i)) / b_alt(i);
        }

        // Update:
        u = u_new;
    }

    // Add u0 at the beginning of u:
    rowvec u_fin(n+1);
    for (int i=0; i<n; i++) {
        u_fin(i+1) = u(i);
    }
    u_fin(0) = u0;
    /* This needs to be done in main.cpp instead of here, or the values of u are
     * changed slightly. BUT it STILL needs to be done here, just without ever
     * using u_fin, so it should not matter at alll, but it does. If else, the
     * values of u are changed slightly...
     * Only God knows why.
     */

    /* The weirdest bug in the universe is happening somewhere in these last few
     * lines.
     */

    return u;
}
