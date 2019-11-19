Imputing missing metabolite values
============================

1. “The original kNN imputation was developed for high-dimensional microarray gene expression data (n «p, n is the number of samples, and p is the number of variables). For each  gene with missing values, this method finds the k nearest genes using Euclidean metric and imputes missing  elements by averaging those non-missing elements of its neighbors. In metabolomics studies, we applied kNN  to find k nearest samples instead and imputed the missing elements. We applied R package impute for this  imputation approach.” (Wei et al 2018, citing Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R, et al. Missing value estimation methods for DNA microarrays. Bioinformatics. 2001;17:520–5.)

2. “QRILC imputation was specifically designed for left-censored data, data missing caused by lower than LOQ. This method imputes missing elements with randomly drawing from a truncated distribution estimated by a quantile regression. A beforehand log-transformation was conducted to improve the imputation accuracy. R package mputeLCMD was applied for this imputation approach. (Wei 2018, citing the impute LCMD package that uses quantile regression).
    - Not clear what exactly is done to perform the regression and imputation. It looks like there are some sophisticated methods in the package, but it doesn’t have a published paper to go along with it.

