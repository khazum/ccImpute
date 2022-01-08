library(SingleCellExperiment)
library(scater)

create_sce_from_counts <- function(counts, colData, rowData = NULL) {
    if(is.null(rowData)) {
        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), 
                                       colData = colData)
    } else {
        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), 
                                       colData = colData,
                                       rowData = rowData)
    }
    # this function writes to logcounts slot
    exprs(sceset) <- log2(calculateCPM(sceset) + 1)
    # use gene names as feature symbols
    rowData(sceset)$feature_symbol <- rownames(sceset)
    # remove features with duplicated names
    if(is.null(rowData)) {
        sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
    }
    # QC
    # isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    # is.spike <- grepl("^ERCC-", rownames(sceset))
    # 
    # sceset <- splitAltExps(sceset, ifelse(is.spike, "ERCC", "gene"))
    # 
    # 
    # perCellQCMetrics(sce,
    #                  exprs_values="counts",
    #                  use_altexps=TRUE)
    # sceset <- calculateQCMetrics(sceset, feature_controls = list("ERCC" = is.spike(sceset, "ERCC")))
    return(sceset)
}

create_sce_from_normcounts <- function(normcounts, colData, rowData = NULL) {
    if(is.null(rowData)) {
        sceset <- SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)), 
                                       colData = colData)
    } else {
        sceset <- SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)), 
                                       colData = colData,
                                       rowData = rowData)
    }
    logcounts(sceset) <- log2(normcounts(sceset) + 1)
    # use gene names as feature symbols
    rowData(sceset)$feature_symbol <- rownames(sceset)
    # remove features with duplicated names
    if(is.null(rowData)) {
        sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
    }
    # QC
    # isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    is.spike <- grepl("^ERCC-", rownames(sceset))
    sce <- splitAltExps(sceset, ifelse(is.spike, "ERCC", "gene"))
    
    return(sceset)
}

create_sce_from_logcounts <- function(logcounts, colData, rowData = NULL) {
    if(is.null(rowData)) {
        sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)), 
                                       colData = colData)
    } else {
        sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)), 
                                       colData = colData,
                                       rowData = rowData)
    }
    # use gene names as feature symbols
    rowData(sceset)$feature_symbol <- rownames(sceset)
    # remove features with duplicated names
    if(is.null(rowData)) {
        sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
    }
    # QC
    # isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    is.spike <- grepl("^ERCC-", rownames(sceset))
    sce <- splitAltExps(sceset, ifelse(is.spike, "ERCC", "gene"))
    return(sceset)
}

