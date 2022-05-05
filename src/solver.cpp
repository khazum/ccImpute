#include <sharedHeader.h>
#include <vector>
#include <unordered_map>

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
#pragma omp parallel for shared(em, col2rows, keys) num_threads(n_cores)
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
            Eigen::VectorXd solution = A.ldlt().solve(b);
            for(unsigned int j = 0; j < row_count; ++j){
                em(rows_indices[j],ci) = solution(j);
            }
        }
    }
    return em;
}
