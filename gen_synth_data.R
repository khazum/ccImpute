suppressPackageStartupMessages({
  library(splatter)
  library(scater)
  library(mclust)
})
source("~/ccImpute/create_sce.R")

for (i in seq(2, 10, by=2)){
  sim <- splatSimulate(group.prob = c(.3, .2, .35, .15), sparsify = FALSE, batchCells=1000*i, nGenes=20000, method = "groups", verbose = FALSE, dropout.type = "experiment")
  cell_type1 <- colData(sim)$Group
  X <- assays(sim)$counts
  sceset <- create_sce_from_counts(X, as.data.frame(cell_type1))
  saveRDS(sceset, paste("~/ccImpute/datasets/sim-", i*1000 , ".rds", sep=""))
}