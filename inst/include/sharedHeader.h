#ifndef cppFunctions_H
#define cppFunctions_H

// [[Rcpp::depends("RcppArmadillo", "RcppEigen")]]
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
#undef Rcpp_hpp
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#endif
