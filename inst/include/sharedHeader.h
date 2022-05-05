#ifndef cppFunctions_H
#define cppFunctions_H

// [[Rcpp::depends("RcppArmadillo", "RcppEigen")]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>

#define NDEBUG 1

#ifdef NDEBUG
#include <RcppEigen.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#endif
