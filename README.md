# ccImpute

## Introduction
Dropout events in scRNA-seq datasets make the lowly expressed genes indistinguishable from true zero expression and different than the low expression present in cells of the same type. This issue makes any subsequent downstream analysis difficult. ccImpute is an imputation algorithm that uses cell similarity established by consensus clustering to impute the most probable dropout events. ccImpute demonstrates performance which exceeds the performance of existing imputation approaches while introducing the least amount of new noise as measured by clustering performance characteristics on datasets with known cell identities. Please visit the [BioConductor page](https://bioconductor.org/packages/release/bioc/html/ccImpute.html) for more info.  

## Installation
To install this package, start R (version "4.2") and enter:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ccImpute")
```

## Issues
Please use [this page](https://github.com/khazum/ccImpute/issues) to report bugs, comments and suggestions.
