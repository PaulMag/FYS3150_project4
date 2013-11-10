//include armadillo;

#ifndef ALGORITHMS_H
#define ALGORITHMSH

arma::rowvec analytic (int, double);

arma::rowvec forwardEuler  (int, double, double);
arma::rowvec backwardEuler (int, double, double);
arma::rowvec crankNicolson (int, double, double);

#endif // ALGORITHMS_H
