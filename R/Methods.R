#' Performs imputation of dropout values in scRNA-seq data using ccImpute
#' algorithm
#'
#' @param logX A normalized and log transformed scRNA-seq expression matrix.
#' @param useRanks A Boolean specifying if non-parametric version of weighted
#' Pearson correlation should be used.
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
#' @param k centers parameter passed to \code{\link[stats]{kmeans}} function
#' @param consMin the low-pass filter threshold for processing consensus
#' matrix.
#' @param nCores the number of cores to be used on the user's machine.
#' If not set, \code{ccImpute} will use all but one cores of your machine.
#' @param kmNStart nstart parameter passed to \code{\link[stats]{kmeans}}.
#' function. Can be set manually. By default it is \code{1000} for up to
#' \code{2000} cells and \code{50} for more than \code{2000} cells.
#' @param kmMax iter.max parameter passed to \code{\link[stats]{kmeans}}.
#'
#' @return A normalized and log transformed scRNA-seq expression matrix with
#' imputed missing values.
#'
#' @importFrom parallel detectCores
#' @importFrom stats prcomp
#' @importFrom matrixStats rowVars
#' @importFrom SIMLR SIMLR_Estimate_Number_of_Clusters
#' @importFrom Rcpp sourceCpp
#' @useDynLib ccImpute
#' @export
impute <- function(logX, useRanks=TRUE, pcaMin, pcaMax, k, consMin=0.65,
                        nCores, kmNStart, kmMax=1000) {

    n <- ncol(logX) #number of samples
    message(c("Running ccImpute on dataset with ", n, " cells."))

    if (missing(nCores)) {# Address missing nCores parameter
        nCores <- parallel::detectCores()
        if (is.null(nCores)) {
            stop("Cannot automatically define a number of available CPU cores.
                    Please specify the nCores parameter manually.")
        }
        if (nCores > 1){nCores <- nCores - 1} # Use all but one core
    }

    # Transform the data to distance matrix form
    distM <- wCorDist(logX, rowVars(logX), useRanks, nCores)
    names(distM) <- c(names, ifelse(useRanks, "Spearman", "Pearson"))

    # Perform a PCA dimensionality reduction on the distance matrix
    distPCA <- prcomp(distM, center = TRUE, scale. = TRUE, retx = FALSE)
    rm(distM) # Conserve space

    # Figure out the sub-datasets in terms of PCA feature indexing
    if(missing(pcaMin) || missing(pcaMax)){
        pcaMin=NULL
        pcaMax=NULL
    }
    nDim <- findNDim(n, distPCA, pcaMin, pcaMax)

    if(missing(k)){  # Address missing k parameter
        k <- SIMLR::SIMLR_Estimate_Number_of_Clusters(logX, NUMC = 2:8)
    }

    # Address missing kmNStart parameter
    if(missing(kmNStart)){kmNStart <- ifelse(n > 2000, 50, 1000)}

    kmResults <- kmeans(distPCA$rotation, k, nCores, nDim, kmNStart, kmMax)
    rm(distPCA) # Conserve space

    consMtx <- getPConsMtx(kmResults, consMin)
    logXt <- t(logX)
    dropIds <- findDropouts(logXt, consMtx)
    impLogX <- t(solveDrops(consMtx, logXt, which(dropIds, arr.ind = TRUE), nCores))
    rownames(impLogX) <- rownames(logX)
    colnames(impLogX) <- colnames(logX)

    message("Imputation finished.")
    return(impLogX)
}

#' This function performs \code{\link[stats]{kmeans}} clustering of the
#' subdatasets corresponding to varying number of PCA loadings as contained in
#' \code{nDim} input parameter.
#'
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
#'
#' @return a list of clustering assignments for all the sub-datasets.
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats kmeans
kmeans <- function(input, k, nCores, nDim, kmNStart, kmMax) {
    cl <- parallel::makeCluster(min(length(nDim), nCores), outfile = "")
    doParallel::registerDoParallel(cl, cores = min(length(nDim), nCores))

    i <- NULL
    results <- foreach::foreach(i = seq_along(nDim), .packages = c("stats")) %dopar% {
        try({
            stats::kmeans(input[, seq_len(nDim[i])], k, iter.max = kmMax,
                            nstart = kmNStart)$cluster
        })
    }

    # stop local cluster
    parallel::stopCluster(cl)

    return(results)
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
    consMtx <- t(apply((consMtx), 2,
                function(i) i/sum(i)))

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
findNDim <- function(n, distPCA, pcaMin, pcaMax){
    isCellCountLow <- n<500
    nDim <- NULL

    # Address missing pcaMin or pcaMax parameters
    if(is.null(pcaMin) || is.null(pcaMax)){
        pcaMin <- ifelse(isCellCountLow, 0.04, 0.0008)
        pcaMax <- ifelse(isCellCountLow, 0.07, 0.002)
    }

    if(isCellCountLow){
        nDim <- floor(pcaMin * n):ceiling(pcaMax * n)
    }
    else{
        eigs <- distPCA$sdev^2
        prop = eigs/sum(eigs)
        n <- length(prop)
        iVal <- n-findInterval(c(pcaMin, pcaMax), prop[n:1])+1
        nDim <- seq(iVal[2], iVal[1], floor((iVal[1]-iVal[2])/15))
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
findDropouts <- function(x, consMtx){
    zeroIndices = x == 0
    xVote <- x
    xVote[zeroIndices] <- -1
    xVote[x > 0] <- 1

    #compute votes - majority wins.
    #Negative value means a true 0 expression, otherwise a dropout.
    votes <- matrix(0L, nrow = nrow(x), ncol = ncol(x))
    votes[zeroIndices] <- (consMtx %*% xVote)[zeroIndices]

    return(votes > 0)
}
