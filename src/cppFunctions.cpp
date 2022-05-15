#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::plugins(openmp)]]

#include <algorithm>
#include <vector>
#include <unordered_map>

using Eigen::Map;
using Eigen::Ref;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;

// Utility function that returns sorted indices rather than elements
VectorXi index_sort(Ref<VectorXd> x, unsigned int s) {
    VectorXi idx(s);
    for(unsigned int i = 0; i < s; ++i){idx[i]=i;}
    std::sort(idx.data(), idx.data()+s, [&](int i, int j){return x[i] < x[j];});
    return idx;
}

// Utility function that returns reverse of sorted indices
VectorXi index_sort_reverse(Ref<VectorXi> x, unsigned int s) {
    VectorXi idx_r(s);
    for(unsigned int i = 0; i < s; ++i){idx_r[x[i]]=i;}
    return idx_r;
}

// Computes weighted ranking of a single vector.
// Helper function needed to compute weighted Spearman correlation.
//
// @param x input vector
// @param w weights for for each element of x
// @return weighted ranking of x
VectorXd wrankFast(Ref<VectorXd> x, Ref<VectorXd> w) {
    unsigned int size = x.size();

    VectorXi ord = index_sort(x, size);
    VectorXi rord = index_sort_reverse(ord, size);

    VectorXd xp = ord.unaryExpr(x);
    VectorXd wp = ord.unaryExpr(w);
    VectorXd rnk(size);
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
    for(unsigned int ii = 0; ii < n; ++ii) {
        rnk[i-ii] = rnki;
    }
    rnk = rord.unaryExpr(rnk);
    return rnk;
}


// Computes weighted ranking of all vectors (columns of x) in parallel
// Helper function needed to compute weighted Spearman correlation.
//
// @param x input vector
// @param w weights for for each element of x
// @return weighted ranking of all vectors (columns of x)
MatrixXd wrankFastAll(Ref<MatrixXd> x, Ref<VectorXd> w, unsigned int n_cores) {
    MatrixXd rank_m(x.rows(), x.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores)
#endif
    for(unsigned int i = 0; i < x.cols(); ++i){
        rank_m.col(i) = wrankFast(x.col(i), w);
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
Eigen::MatrixXd wCorDist(const Eigen::Map<Eigen::MatrixXd> x, 
                         const Eigen::Map<Eigen::VectorXd> w,
                         const bool useRanks,
                   const unsigned int n_cores){
    
    Eigen::MatrixXd temp = useRanks ? wrankFastAll(x, w, n_cores) : x;
    double sumw = w.sum();
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores) 
    #endif
    for(unsigned int i = 0; i < temp.cols(); ++i){
        Eigen::ArrayXXd col_arr = temp.col(i).array();
        temp.col(i) = col_arr - (w.array() * col_arr).sum()/sumw;
    }
    
    Eigen::MatrixXd corrs(temp.cols(),temp.cols());
    Eigen::ArrayXXd w_arr = w.array();
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_cores) schedule(dynamic)
    #endif
    for(unsigned int i = 0; i < temp.cols(); ++i){
        Eigen::ArrayXXd t1 = temp.col(i).array();
        for(unsigned int j = i+1; j < temp.cols(); ++j){
            Eigen::ArrayXXd t2 = temp.col(j).array();
            double n = (w_arr*t1*t2).sum();
            double d = (w_arr*t1.pow(2)).sum() * (w_arr*t2.pow(2)).sum();
            // double d = p1*p2;
            d = pow(d, 0.5);
            corrs(i,j) = 1 - n/d;
            corrs(j,i) = corrs(i,j);
        }
    }
    return(corrs);
}

//' Computes consensus matrix given cluster labels
//'
//' @param dat a matrix containing clustering solutions in columns
//' @return consensus matrix
// [[Rcpp::export]]
Eigen::MatrixXd getConsMtx(const Eigen::Map<Eigen::MatrixXi> dat) {
    Eigen::MatrixXd res = dat.cols() * 
        Eigen::MatrixXd::Identity(dat.rows(), dat.rows());
    unsigned int i, j, k;
    for (unsigned j = 0; j < dat.cols(); ++j) {
        for (unsigned i = 0; i < dat.rows(); ++i) {
            for (unsigned k = i + 1; k < dat.rows(); ++k) {
                if (dat(i, j) == dat(k, j)) {
                    ++res(i, k);
                    ++res(k, i);
                }
            }
        }
    }
    res /= dat.cols();
    return res;
}

//' Computes imputed expression matrix using linear eq solver.
//'
//'@param cm processed consensus matrix
//'@param em expression matrix
//'@param ids location of values determined to be dropout events
//'@param n_cores number of cores to use for parallel computation.
//'@return imputed expression matrix
// [[Rcpp::export]]
Eigen::MatrixXd solveDrops(const Eigen::Map<Eigen::MatrixXd> cm, 
                           Eigen::Map<Eigen::MatrixXd> em, 
                           const Eigen::Map<Eigen::MatrixXi> ids, 
                           const int n_cores) {
    //map row to columns
    std::unordered_map<unsigned int, std::vector<unsigned int>> col2rows;
    std::vector<unsigned int> keys;
    
    // Map out the rows for each column
    for(unsigned int i=0; i<ids.rows(); ++i){
        // account for the fact that r indices are 1 based
        unsigned int col = ids(i,1) - 1;
        if (col2rows.find(col) == col2rows.end()) {
            keys.push_back(col);
            col2rows[col] = std::vector<unsigned int>();
        }
        // account for the fact that r indices are 1 based
        col2rows[col].push_back(ids(i,0)-1);
    }
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores)
#endif
    for(unsigned int i=0; i < keys.size(); ++i){
        unsigned int ci = keys[i];
        std::vector<unsigned int> rows_indices = col2rows[ci];
        unsigned int row_count = rows_indices.size();
        if(row_count == 1){
            em(rows_indices[0],ci) =cm.row(rows_indices[0]).dot(em.col(ci));
        }
        else{
            Eigen::MatrixXd A(row_count, row_count);
            
            for(unsigned int j=0; j < row_count; ++j){
                for(unsigned int k=0; k < row_count; ++k){
                    A(j,k) = j==k ? 1:-cm(rows_indices[j], rows_indices[k]);
                }
            }
            Eigen::VectorXd b(row_count);
            for(unsigned int j=0; j < row_count; ++j){
                b[j] = cm.row(rows_indices[j]).dot(em.col(ci));
            }
            Eigen::VectorXd solution = A.llt().solve(b);
            for(unsigned int j = 0; j < row_count; ++j){
                em(rows_indices[j],ci) = solution(j);
            }
        }
    }
    return em;
}