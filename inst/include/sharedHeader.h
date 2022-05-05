#ifndef cppFunctions_H
#define cppFunctions_H

// [[Rcpp::depends("RcppArmadillo", "RcppEigen")]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>

#define NDEBUG 
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#endif