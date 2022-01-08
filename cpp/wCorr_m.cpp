#define ARMA_ALLOW_FAKE_GCC
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <Rmath.h>
#include <omp.h>

using namespace Rcpp;

arma::vec wrankFast(arma::vec x, const arma::vec& w) {
  const arma::uvec ord = arma::sort_index(x);
  int size = ord.size();
  arma::vec temp = arma::linspace(0, size-1, size);
  const arma::uvec ind = arma::sort_index(ord);
  
  temp = temp.elem(ind);
  arma::uvec rord = arma::conv_to<arma::uvec>::from(temp);
  const arma::vec xp = x.elem(ord);
  const arma::vec wp = w.elem(ord);
  arma::vec rnk(size);
  double t1 = 0;
  int i = 0;
  double t2 = 0;
  double n = 0 ;
  while(i < size-1) {
    
    t2 = t2 + wp[i];
    n = n + 1;
    if(xp[i+1] != xp[i]) {
      double rnki = t1 + (n+1)/(2*n)*t2;
      for(int ii =0; ii < n; ii++) {
        rnk[i-ii] = rnki;
        
      }
      t1 = t1 + t2; 
      t2 = 0;
      n = 0;
    }
    i = i + 1;
  }
  n = n + 1;
  t2 = t2 + wp[i];
  double rnki = t1 + (n+1)/(2*n)*t2; 
  for(int ii = 0; ii < n; ii++) {
    rnk[i-ii] = rnki;
  }
  rnk = rnk.elem(rord);
  return(rnk);
}

double cont(arma::vec x,  arma::vec y, const arma::vec& w) {
  double sumw = arma::sum(w);
  double xb = arma::sum(w%x)/sumw;
  double yb = arma::sum(w%y)/sumw;
  const arma::vec temp1 = x-xb;
  const arma::vec temp2 = y-yb;
  double numerator = arma::sum(w%temp1%temp2);
  double denominator = pow(arma::sum(w%arma::square(temp1)) * sum(w%arma::square(temp2)), 0.5);
  return numerator/denominator;
}   

arma::mat wrankFast_m(arma::mat x, const arma::vec& w) {
  arma::mat rank_m(x.n_rows, x.n_cols);
  int procs = omp_get_num_procs();
  omp_set_num_threads(procs-1);
  #pragma omp parallel for shared(rank_m)
  for(int i = 0; i < x.n_cols; ++i){
    rank_m.col(i) = wrankFast(x.col(i),w);
  }
  return(rank_m);
}

// [[Rcpp::export]]
arma::mat w_cor_dist_m(arma::mat x, const arma::vec& w){
  arma::mat rank_m(x.n_rows, x.n_cols);
  int procs = omp_get_num_procs();
  omp_set_num_threads(procs-1);

  arma::mat temp = wrankFast_m(x,w);
  
  double sumw = arma::sum(w);
  
  #pragma omp parallel for shared(temp)
  for(int i = 0; i < temp.n_cols; ++i){
    temp.col(i) = temp.col(i) - arma::sum(w%temp.col(i))/sumw;
  }
  
  arma::mat correlations(temp.n_cols,temp.n_cols);
  
  #pragma omp parallel for shared(correlations)
  for(int i = 0; i < temp.n_cols; ++i){
    for(int j = i+1; j < temp.n_cols; ++j){
      arma::vec temp1 = temp.col(i);
      arma::vec temp2 = temp.col(j);
      double numerator = arma::sum(w%temp1%temp2);
      double denominator = pow(arma::sum(w%arma::square(temp1)) * sum(w%arma::square(temp2)), 0.5);
      correlations(i,j) = 1 - numerator/denominator;
      correlations(j,i) = correlations(i,j);
    }
  }
  return(correlations);
}