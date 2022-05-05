#include <RcppArmadillo.h>
#define NDEBUG 
#include <RcppEigen.h>

#include <Rmath.h>
#include <cstdlib>
#include <Rmath.h>
#include <vector>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends("RcppArmadillo", "RcppEigen")]]
// [[Rcpp::plugins(openmp)]]


// Computes weighted ranking of a single vector.
// Helper function needed to compute weighted Spearman correlation.
//
// @param x input vector
// @param w weights for for each element of x
// @return weighted ranking of x
arma::vec wrankFast(const arma::vec& x, const arma::vec& w) {
    const arma::uvec ord = arma::sort_index(x);
    unsigned int size = ord.size();
    arma::vec temp = arma::linspace(0, size-1, size);
    const arma::uvec ind = arma::sort_index(ord);
    
    temp = temp.elem(ind);
    arma::uvec rord = arma::conv_to<arma::uvec>::from(temp);
    const arma::vec xp = x.elem(ord);
    const arma::vec wp = w.elem(ord);
    arma::vec rnk(size);
    double t1 = 0;
    unsigned int i = 0;
    double t2 = 0;
    double n = 0 ;
    while(i < size-1) {
        t2 = t2 + wp[i];
        n = n + 1;
        if(xp[i+1] != xp[i]) {
            double rnki = t1 + (n+1)/(2*n)*t2;
            for(int ii =0; ii < n; ++ii) {
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
    for(unsigned int ii = 0; ii < n; ii++) {
        rnk[i-ii] = rnki;
    }
    rnk = rnk.elem(rord);
    return rnk;
}

// Computes weighted ranking of all vectors (columns of x) in parallel
// Helper function needed to compute weighted Spearman correlation.
//
// @param x input vector
// @param w weights for for each element of x
// @return weighted ranking of all vectors (columns of x)
arma::mat wrankFastAll(const arma::mat& x, const arma::vec& w, 
                           unsigned int n_cores) {
    arma::mat rank_m(x.n_rows, x.n_cols);
#ifdef _OPENMP
#pragma omp parallel for shared(rank_m)  num_threads(n_cores)
#endif
    for(unsigned int i = 0; i < x.n_cols; ++i){
        rank_m.col(i) = wrankFast(x.col(i),w);
    }
    return rank_m;
}

//' Computes a weighted Pearson distance measure matrix. If ranks are used
//' this measure turns into weighted Spearman distance measure matrix.
//'
//'@param x input with columns containing each observation
//'@param w weights for all values in a obervation
//'@param useRanks indicates if Pearson should be computed on weighted ranks.
//'@param n_cores number of cores to use for parallel computation.
//'@return weighted Pearson distance measure matrix. If ranks are used
//' this measure turns into weighted Spearman distance measure matrix.
// [[Rcpp::export]]
arma::mat wCorDist(arma::mat x, const arma::vec& w, const bool useRanks,
                       const unsigned int n_cores){
    arma::mat rank_m(x.n_rows, x.n_cols);
    
    arma::mat temp = useRanks ? wrankFastAll(x, w, n_cores) : x;
    double sumw = arma::sum(w);
    
#ifdef _OPENMP
#pragma omp parallel for shared(temp) num_threads(n_cores)
#endif
    for(unsigned int i = 0; i < temp.n_cols; ++i){
        temp.col(i) = temp.col(i) - arma::sum(w%temp.col(i))/sumw;
    }
    arma::mat corrs(temp.n_cols,temp.n_cols);
    
#ifdef _OPENMP
#pragma omp parallel for shared(corrs) num_threads(n_cores)
#endif
    for(unsigned int i = 0; i < temp.n_cols; ++i){
        for(unsigned int j = i+1; j < temp.n_cols; ++j){
            arma::vec t1 = temp.col(i);
            arma::vec t2 = temp.col(j);
            double n = arma::sum(w%t1%t2);
            double d = arma::sum(w%arma::square(t1)) * sum(w%arma::square(t2));
            d = pow(d, 0.5);
            corrs(i,j) = 1 - n/d;
            corrs(j,i) = corrs(i,j);
        }
    }
    return(corrs);
}