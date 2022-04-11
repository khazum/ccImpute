#define ARMA_ALLOW_FAKE_GCC
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <omp.h>
#include <cstdlib>
#include <vector>
#include <unordered_map>

using namespace arma;

// [[Rcpp::plugins(openmp)]]

//' Computes consensus matrix given cluster labels
//'
//' @param dat a matrix containing clustering solutions in columns
//' @return consensus matrix
// [[Rcpp::export]]
arma::mat getConsMtx(const arma::mat dat) {

    mat res = dat.n_cols * eye<mat>(dat.n_rows, dat.n_rows);

    unsigned int i, j, k;
    for (j = 0; j < dat.n_cols; j++) {
        for (i = 0; i < dat.n_rows; i++) {
            for (k = i + 1; k < dat.n_rows; k++) {
                if (dat(i, j) == dat(k, j)) {
                    ++res(i, k);
                    ++res(k, i);
                }
            }
        }
    }
    res /= dat.n_cols;
    return res;
}

// Computes weighted ranking of a single vector.
// Helper function needed to compute weighted Spearman correlation.
//
// @param x input vector
// @param w weights for for each element of x
// @return weighted ranking of x
arma::vec wrankFast(arma::vec x, const arma::vec& w) {
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
    return(rnk);
}

// Computes weighted ranking of all vectors (columns of x) in parallel
// Helper function needed to compute weighted Spearman correlation.
//
// @param x input vector
// @param w weights for for each element of x
// @return weighted ranking of all vectors (columns of x)
arma::mat wrankFastAll(arma::mat x, const arma::vec& w) {
    arma::mat rank_m(x.n_rows, x.n_cols);
    #pragma omp parallel for shared(rank_m)
    for(unsigned int i = 0; i < x.n_cols; ++i){
        rank_m.col(i) = wrankFast(x.col(i),w);
    }
    return(rank_m);
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
    omp_set_num_threads(n_cores);

    arma::mat temp = useRanks ? wrankFastAll(x,w) : x;

    double sumw = arma::sum(w);

    #pragma omp parallel for shared(temp)
    for(unsigned int i = 0; i < temp.n_cols; ++i){
        temp.col(i) = temp.col(i) - arma::sum(w%temp.col(i))/sumw;
    }
    arma::mat corrs(temp.n_cols,temp.n_cols);

    #pragma omp parallel for shared(corrs)
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

//' Computes imputed expression matrix using linear eq solver.
//'
//'@param cm processed consensus matrix
//'@param em expression matrix
//'@param ids location of values determined to be dropout events
//'@param n_cores number of cores to use for parallel computation.
//'@return imputed expression matrix
// [[Rcpp::export]]
arma::mat solveDrops(arma::mat cm, arma::mat em, arma::mat ids, const int n_cores) {
    omp_set_num_threads(n_cores);

    //map row to columns
    std::unordered_map<int, std::vector<int>> col2rows;
    std::vector<int> keys;

    // Map out the rows for each column
    for(unsigned int i=0; i<ids.n_rows; ++i){
        // account for the fact that r indices are 1 based
        int col = ids.at(i,1) - 1;
        if (col2rows.find(col) == col2rows.end()) {
            keys.push_back(col);
            col2rows[col] = std::vector<int>();
        }
        // account for the fact that r indices are 1 based
        col2rows[col].push_back(ids.at(i,0)-1);
    }

    #pragma omp parallel for shared(em, col2rows, keys)
    for(unsigned int i=0; i < keys.size(); ++i){
        unsigned int ci = keys[i];
        std::vector<int> rows_indices = col2rows[ci];
        unsigned int row_count = rows_indices.size();
        if(row_count == 1){
            em.at(rows_indices[0],ci) = dot(cm.row(rows_indices[0]), em.col(ci));
        }
        else{
            arma::mat A(row_count, row_count);

            for(unsigned int j=0; j < row_count; ++j){
                for(unsigned int k=0; k < row_count; ++k){
                    A.at(j,k) = j==k ? 1:-cm.at(rows_indices[j], rows_indices[k]);
                }
            }

            arma::vec b(row_count);
            for(unsigned int j=0; j < row_count; ++j){
                b[j] = dot(cm.row(rows_indices[j]), em.col(ci));
            }
            arma::vec solution = solve(A, b, arma::solve_opts::fast);

            for(unsigned int j = 0; j < row_count; ++j){
                em.at(rows_indices[j],ci) = solution.at(j);
            }
        }
    }
    return em;
}
