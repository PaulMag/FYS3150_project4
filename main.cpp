#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

#include "algorithms.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {

    istringstream nS  (argv[1]);
    istringstream dtFactorS (argv[2]);

    int n;
    double dtFactor;
    nS        >> n;
    dtFactorS >> dtFactor;

    double dt = dtFactor / (n*n); // makes sure forwardEuler() is stable
    double time1 = 0.05;
    double time2 = 1.0;

    /* I want to compare each of the 3 solving methods with the analytical
     * solution for a short time and for a time long enough for the system to
     * reach its stationary state.
     */

    rowvec anal1 =      analytic(n, time1);
    rowvec anal2 =      analytic(n, time2);

    rowvec forw1 =  forwardEuler(n, time1, dt);
    rowvec forw2 =  forwardEuler(n, time2, dt);

    rowvec bacw1 = backwardEuler(n, time1, dt);
    rowvec bacw2 = backwardEuler(n, time2, dt);

    rowvec crni1 = crankNicolson(n, time1, dt);
    rowvec crni2 = crankNicolson(n, time2, dt);

    // begin bugfix
    /* Add 1.0 at the beginning of the Crank-Nicolson methods so the length of
     * these vectors becomes n+1 instead of n. This should of course be done
     * inside the function and not here, but that didn't work, because of the
     * weirdest bug in the universe, AKA I have no idea why.
     */
    rowvec crni1_fin(n+1);
    rowvec crni2_fin(n+1);
    for (int i=0; i<n; i++) {
        crni1_fin(i+1) = crni1(i);
        crni2_fin(i+1) = crni2(i);
    }
    crni1_fin(0) = 1.0;
    crni2_fin(0) = 1.0;
    // end bugfix

    for (int i=0; i<n; i++) {
        forw1(i) = abs(forw1(i) - anal1(i));
        forw2(i) = abs(forw2(i) - anal2(i));
        bacw1(i) = abs(bacw1(i) - anal1(i));
        bacw2(i) = abs(bacw2(i) - anal2(i));
        crni1_fin(i) = abs(crni1_fin(i) - anal1(i));
        crni2_fin(i) = abs(crni2_fin(i) - anal2(i));
    }

    ofstream outfile;
    outfile.open("data/error.dat");
    outfile << n << " " << time1 << " " << time2 << " " << dt << endl;
    outfile << forw1 << forw2 << bacw1 << bacw2 << crni1_fin << crni2_fin;
    outfile.close();

}
