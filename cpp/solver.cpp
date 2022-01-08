#define ARMA_ALLOW_FAKE_GCC
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <omp.h>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <omp.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat solve_dropouts(arma::mat cons_m, arma::mat exp_m, arma::mat indices) {
  //map row to columns
  std::unordered_map<int, std::vector<int>> col2rows;
  std::vector<int> keys;
  
  // Map out the rows for each column
  for(int i=0; i<indices.n_rows; ++i){
    int col = indices.at(i,1) - 1; // account for the fact that r indices are 1 based
    if (col2rows.find(col) == col2rows.end()) {
      keys.push_back(col);
      col2rows[col] = std::vector<int>();
    }
    col2rows[col].push_back(indices.at(i,0)-1); // account for the fact that r indices are 1 based
  }
  
#pragma omp parallel for shared(exp_m, col2rows, keys)
  for(int i=0; i < keys.size(); ++i){
    int col_index = keys[i];
    std::vector<int> rows_indices = col2rows[col_index];
    int row_count = rows_indices.size();
    if(row_count == 1){
      exp_m.at(rows_indices[0],col_index) = dot(cons_m.row(rows_indices[0]), exp_m.col(col_index));
    }
    else{
      arma::mat A(row_count, row_count);
      
      for(int j=0; j < row_count; ++j){
        for(int k=0; k < row_count; ++k){
          A.at(j,k) = j==k ? 1:-cons_m.at(rows_indices[j], rows_indices[k]);
        }
      }
      
      arma::vec b(row_count);
      for(int j=0; j < row_count; ++j){
        b[j] = dot(cons_m.row(rows_indices[j]), exp_m.col(col_index));
      }
      arma::vec solution = solve(A, b, arma::solve_opts::fast);
      
      for(int j = 0; j < row_count; ++j){
        exp_m.at(rows_indices[j],col_index) = solution.at(j);
      }
    }
  }
  return exp_m;
}