#' Performs imputation of dropout values in scRNA-seq data using ccImpute
#' algorithm as described in the ccImpute: an accurate and scalable consensus 
#' clustering based algorithm to impute dropout events in the single-cell 
#' RNA-seq data DOI: https://doi.org/10.1186/s12859-022-04814-8
#'
#' @param logX A normalized and log transformed scRNA-seq expression matrix.
#' @param useRanks A Boolean specifying if non-parametric version of weighted
#' Pearson correlation should be used. It's recommended to keep this as TRUE
#' since  this performs better as determined experimentally. However, FALSE
#' will also provide decent results with the benefit or faster runtime. 
#' @param pcaMin This is used to establish the number of minimum PCA features
#' used for generating subsets. For small datasets up to \code{500} cells this
#' equals pcaMin*n minimum features, where \code{n} is number of cells. For
#' large datasets, this corresponds to the feature count that has proportion of
#' variance less than \code{pcaMin}. Both pcaMin and pcaMax must be specified
#' to be considered. It's best to keep this value as default unless a better
#' value was obtained experimentally.
#' @param pcaMax This is used to establish the number of maximum PCA features
#' used for generating subsets. For small datasets up to \code{500} cells this
#' equals pcaMax*n maximum features, where \code{n} is number of cells. For
#' large datasets, this corresponds to the feature count that has proportion of
#' variance less than \code{pcaMax}. Both pcaMin and pcaMax must be specified
#' to be considered. It's best to keep this value as default unless a better
#' value was obtained experimentally.
#' @param k centers parameter passed to \code{\link[stats]{kmeans}} function.
#' This corresponds to a number of different cell groups in data. This can be
#' estimated in a number of methods. If not provided we take the approach
#' provided in the SIMLR package. 
#' (https://www.bioconductor.org/packages/release/bioc/html/SIMLR.html)
#' @param consMin the low-pass filter threshold for processing consensus
#' matrix. This is to eliminate noise from unlikely clustering assignmnets.
#' It is recommended to keep this value >-.5.
#' @param nCores the number of cores to be used on the user's machine.
#' If not set, \code{ccImpute} will use all but one cores of your machine.
#' @param kmNStart nstart parameter passed to \code{\link[stats]{kmeans}}.
#' function. Can be set manually. By default it is \code{1000} for up to
#' \code{2000} cells and \code{50} for more than \code{2000} cells.
#' @param kmMax iter.max parameter passed to \code{\link[stats]{kmeans}}.
#' @param rand_seed sets the seed of the random number generator.
#' \code{ccImpute} is a stochastic method, and setting the \code{rand_seed}
#' allows reproducibility.
#'
#' @return A normalized and log transformed scRNA-seq expression matrix with
#' imputed missing values.
#'
#' @importFrom stats prcomp
#' @importFrom matrixStats rowVars
#' @importFrom Rcpp sourceCpp
#' @importFrom SIMLR SIMLR_Estimate_Number_of_Clusters
#' @importFrom BiocParallel bplapply bpnworkers bpparam
#'
#' @useDynLib ccImpute
#' @export
#' @examples
#' exp_matrix <- log(abs(matrix(rnorm(1000000),nrow=10000))+1)
#' ccImpute(exp_matrix, k = 2)
ccImpute <- function(logX, useRanks=TRUE, pcaMin, pcaMax, k, consMin=0.65,
                        kmNStart, kmMax=1000) {
    logX <- as.matrix(logX) 
    n <- ncol(logX) #number of samples
    message(c("Running ccImpute on dataset with ", n, " cells."))

    nCores <- bpnworkers(bpparam())

    if(missing(k)){
        k <- SIMLR::SIMLR_Estimate_Number_of_Clusters(logX, NUMC = 2:12)
        k <- which.min((k$K1+k$K2))+1
    }

    # Transform the data to distance matrix form
    distM <- wCorDist(logX, rowVars(logX), useRanks, nCores)
    names(distM) <- c(names, ifelse(useRanks, "Spearman", "Pearson"))

    # Perform a PCA dimensionality reduction on the distance matrix
    distPCA <- prcomp(distM, center = TRUE, scale. = TRUE, retx = FALSE)
    rm(distM) # Conserve space

    # Figure out the sub-datasets in terms of PCA feature indexing
    if(missing(pcaMin) || missing(pcaMax)){pcaMin <- pcaMax <- NULL}
    nDim <- findNDim(n, distPCA, pcaMin, pcaMax)

    # Address missing kmNStart parameter
    if(missing(kmNStart)){kmNStart <- ifelse(n > 2000, 50, 1000)}

    kmResults <- bplapply(nDim, kmAux, distPCA$rotation[], k, kmNStart, kmMax)
    rm(distPCA) # Conserve space

    consMtx <- getPConsMtx(kmResults, consMin)
    logXt <- t(logX)
    dropIds <- findDropouts(logXt, consMtx)
    iLogX <- t(solveDrops(consMtx, logXt, which(dropIds, arr.ind = TRUE),
                            nCores))

    rownames(iLogX) <- rownames(logX)
    colnames(iLogX) <- colnames(logX)
    message("Imputation finished.")
    return(iLogX)
}

#' This function performs \code{\link[stats]{kmeans}} clustering of the
#' subdataset corresponding to a given range i of PCA loadings as contained in
#' \code{input} parameter.
#' 
#' @param i number of loadings to use.
#' @param input the matrix of all variable loadings.
#' @param k centers (integer) parameter passed to \code{\link[stats]{kmeans}}
#' function.
#' @param nCores defines the number of cores to be used on the user's machine.
#' If not set, `ccImpute` will use all but one cores of your machine.
#' @param nDim the list of containing a number of PCA loadings to use for each
#' sub-dataset.
#' @param kmNStart nstart parameter passed to \code{\link[stats]{kmeans}}
#' function. Can be set manually. By default it is \code{1000} for up to
#' \code{2000} cells and \code{50} for more than \code{2000} cells.
#' @param kmMax iter.max parameter passed to
#' \code{\link[stats]{kmeans}} function.
#' @param rand_seed sets the seed of the random number generator.
#' \code{ccImpute} is a stochastic method, and setting the \code{rand_seed}
#' allows reproducibility.
#'
#' @return a list of clustering assignments for all the sub-datasets.
#' @keywords internal
#' @importFrom stats kmeans
kmAux <- function(i, input, k, kmNStart, kmMax) {
    x <- input[, seq_len(i)]
    return(stats::kmeans(x, k, iter.max = kmMax, nstart = kmNStart)$cluster)
}

#' Get processed consensus matrix.
#'
#' This function gets consensus matrix based on the clustering solutions
#' contained in the \code{kmResults} input parameter and does the processing to
#' use it for imputation.
#'
#' @param kmResults list of k-means clustering assignments on the PCA loadings
#' sub-datasets.
#' @param consMin the low-pass filter threshold value for processed
#' consensus matrix.
#'
#' @return a processed consensus matrix.
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @useDynLib ccImpute
getPConsMtx <- function(kmResults, consMin) {
    d <- matrix(unlist(kmResults), nrow = length(kmResults[[1]]))
    consMtx <- getConsMtx(d)

    # Remove diagonal entries
    consMtx <- consMtx - diag(nrow(consMtx))

    # Low pass filter
    consMtx[consMtx < consMin] <- 0

    # Normalize the entries
    consMtx <- t(apply((consMtx), 2, function(i) i/sum(i)))

    # Replace NA values with 0
    consMtx[is.na(consMtx)] <- 0
    return(consMtx)
}

#' Establish what subsets of loadings from PCA distance measure are used for
#' for measuring cluster instability
#'
#' @param n number of samples
#' @param distPCA PCA reduced distance matrix
#' @param pcaMin This is used to establish the number of minimum PCA features
#' used for generating subsets. For small datasets up to \code{500} cells this
#' equals pcaMin*n minimum features, where \code{n} is number of cells. For
#' large datasets, this corresponds to the feature count that has proportion of
#' variance less than \code{pcaMin}. Both pcaMin and pcaMax must be specified
#' to be considered.
#' @param pcaMax This is used to establish the number of maximum PCA features
#' used for generating subsets. For small datasets up to \code{500} cells this
#' equals pcaMax*n maximum features, where \code{n} is number of cells. For
#' large datasets, this corresponds to the feature count that has proportion of
#' variance less than \code{pcaMax}. Both pcaMin and pcaMax must be specified
#' to be considered.
#' @return list of numbers with each number corresponding to the number of
#' loadings to use for clustering.
#' @keywords internal
findNDim <- function(n, distPCA, pcaMin, pcaMax){
    isCellCountLow <- n<500
    nDim <- NULL

    # Address missing pcaMin or pcaMax parameters
    if(is.null(pcaMin) || is.null(pcaMax)){
        pcaMin <- ifelse(isCellCountLow, 0.04, 0.0008)
        pcaMax <- ifelse(isCellCountLow, 0.07, 0.002)
    }

    if(isCellCountLow){
        lv <- floor(pcaMin * n)
        rv <- ceiling(pcaMax * n)
        nDim <- seq(from=lv,to=rv, by=ceiling((rv-lv)/14))
    }
    else{
        eigs <- distPCA$sdev^2
        prop <- eigs/sum(eigs)
        iVal <- n-findInterval(c(pcaMin, pcaMax), prop[n:1])+1
        nDim <- seq(iVal[2], iVal[1], ceiling((iVal[1]-iVal[2])/14))
    }
    return(nDim)
}

#' Establishes which zero values in \code{x} are dropout events based on
#' weighted cell voting with weights derived from processed consensus matrix
#' \code{consMtx}.
#'
#' @param x transpose of log normalized expression matrix
#' @param consMtx processed consensus matrix
#' 
#' @return list of indices in x that are dropout events
#' @keywords internal
findDropouts <- function(x, consMtx){
    zeroIndices <- x == 0
    xVote <- x
    xVote[zeroIndices] <- -1
    xVote[x > 0] <- 1

    #compute votes - majority wins.
    #Negative value means a true 0 expression, otherwise a dropout.
    votes <- matrix(0L, nrow = nrow(x), ncol = ncol(x))
    votes[zeroIndices] <- (consMtx %*% xVote)[zeroIndices]

    return(votes > 0)
}
