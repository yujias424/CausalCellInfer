# CausalCellInfer

### Description

This is the package used for deconvoluting the bulk RNA-seq data using the single-cell RNA-seq data as reference. The idea is inherited from the invariant causal prediction (ICP). 
We utilized the ICP approach to identify the stable markers for each cells and using the deep-learning approach to deconvolute the bulk RNA-seq data. Additionally, we add several 
analysis module in this package, including inferring the cell-type specific expression profile (CTS) using ENIGMA and identifying the causal genes using our recent proposed method.

### Installation

Please use the devtools::install_github to install the package.
```
library(devtools)
install_github("yujias424/CausalCellInfer")
```

### Usage

Please refer to the vignettes.

### Required Package

In brief, you need to install the following packages: 

1. R packages:

    ENIGMA, abind, data.table, pcalg, reticulate, scater, WGCNA, coop, reticulate, stats, xgboost, zellkonverter
2. Python packages:

    torch, numpy, pandas, matplotlib, seaborn, anndata, sklearn

We strongly recommend you to install the pcalg we provided in the Github repository, as we have revise some of the functions provided in the orignal package.

Additionally, we attached the output of sessionInfo() to show the exact version of the package we used in our analysis.

```
> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: x86_64-conda-linux-gnu
Running under: Ubuntu 24.04 LTS

Matrix products: default
BLAS/LAPACK: /home/yujia/miniconda3/envs/causalcellinfer/lib/libopenblasp-r0.3.21.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Hong_Kong
tzcode source: system (glibc)

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] preprocessCore_1.66.0       purrr_1.0.2                
 [3] tibble_3.2.1                nnls_1.5                   
 [5] sva_3.52.0                  BiocParallel_1.38.0        
 [7] genefilter_1.86.0           mgcv_1.9-1                 
 [9] nlme_3.1-165                doParallel_1.0.17          
[11] iterators_1.0.14            foreach_1.5.2              
[13] CausalCellInfer_1.0         scater_1.32.1              
[15] ggplot2_3.5.1               scuttle_1.14.0             
[17] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
[19] Biobase_2.64.0              GenomicRanges_1.56.1       
[21] GenomeInfoDb_1.40.1         IRanges_2.38.1             
[23] S4Vectors_0.42.1            BiocGenerics_0.50.0        
[25] MatrixGenerics_1.16.0       matrixStats_1.4.1          
[27] zellkonverter_1.14.1        stringr_1.5.1              
[29] xgboost_1.7.8.1             data.table_1.15.4          

loaded via a namespace (and not attached):
  [1] splines_4.4.1             filelock_1.0.3           
  [3] R.oo_1.26.0               basilisk.utils_1.16.0    
  [5] graph_1.82.0              XML_3.99-0.17            
  [7] rpart_4.1.23              lifecycle_1.0.4          
  [9] fastcluster_1.2.6         edgeR_4.2.1              
 [11] lattice_0.22-6            MASS_7.3-61              
 [13] backports_1.5.0           magrittr_2.0.3           
 [15] limma_3.60.6              Hmisc_5.1-3              
 [17] rmarkdown_2.28            reticulate_1.39.0        
 [19] DBI_1.2.3                 abind_1.4-5              
 [21] zlibbioc_1.50.0           sfsmisc_1.1-19           
 [23] R.utils_2.12.3            nnet_7.3-19              
 [25] GenomeInfoDbData_1.2.12   ggrepel_0.9.6            
 [27] irlba_2.3.5.1             ENIGMA_0.1.6             
 [29] annotate_1.82.0           DelayedMatrixStats_1.26.0
 [31] codetools_0.2-20          DelayedArray_0.30.1      
 [33] tidyselect_1.2.1          UCSC.utils_1.0.0         
 [35] ScaledMatrix_1.12.0       viridis_0.6.5            
 [37] dynamicTreeCut_1.63-1     base64enc_0.1-3          
 [39] jsonlite_1.8.9            BiocNeighbors_1.22.0     
 [41] e1071_1.7-16              Formula_1.2-5            
 [43] survival_3.7-0            tools_4.4.1              
 [45] Rcpp_1.0.13               glue_1.8.0               
 [47] gridExtra_2.3             SparseArray_1.4.8        
 [49] xfun_0.48                 HDF5Array_1.32.1         
 [51] dplyr_1.1.4               withr_3.0.1              
 [53] fastmap_1.2.0             basilisk_1.16.0          
 [55] rhdf5filters_1.16.0       fansi_1.0.6              
 [57] digest_0.6.37             rsvd_1.0.5               
 [59] R6_2.5.1                  colorspace_2.1-1         
 [61] GO.db_3.19.1              RSQLite_2.3.7            
 [63] R.methodsS3_1.8.2         utf8_1.2.4               
 [65] tidyr_1.3.1               generics_0.1.3           
 [67] pcalg_2.6-12              corpcor_1.6.10           
 [69] robustbase_0.99-4-1       class_7.3-22             
 [71] httr_1.4.7                htmlwidgets_1.6.4        
 [73] S4Arrays_1.4.1            pkgconfig_2.0.3          
 [75] gtable_0.3.5              blob_1.2.4               
 [77] impute_1.78.0             XVector_0.44.0           
 [79] htmltools_0.5.8.1         RBGL_1.80.0              
 [81] clue_0.3-65               scales_1.3.0             
 [83] png_0.1-8                 knitr_1.48               
 [85] rstudioapi_0.16.0         ggm_2.3                  
 [87] reshape2_1.4.4            checkmate_2.3.2          
 [89] bdsmatrix_1.3-7           rhdf5_2.48.0             
 [91] proxy_0.4-27              cachem_1.1.0             
 [93] vipor_0.4.7               foreign_0.8-87           
 [95] AnnotationDbi_1.66.0      pillar_1.9.0             
 [97] grid_4.4.1                fastICA_1.2-5.1          
 [99] vctrs_0.6.5               BiocSingular_1.20.0      
[101] beachmat_2.20.0           xtable_1.8-4             
[103] cluster_2.1.6             beeswarm_0.4.0           
[105] htmlTable_2.4.3           evaluate_1.0.0           
[107] cli_3.6.3                 locfit_1.5-9.10          
[109] compiler_4.4.1            rlang_1.1.4              
[111] crayon_1.5.3              plyr_1.8.9               
[113] ggbeeswarm_0.7.2          stringi_1.8.4            
[115] viridisLite_0.4.2         WGCNA_1.73               
[117] munsell_0.5.1             Biostrings_2.72.1        
[119] coop_0.6-3                Matrix_1.7-0             
[121] dir.expiry_1.12.0         sparseMatrixStats_1.16.0 
[123] bit64_4.5.2               Rhdf5lib_1.26.0          
[125] KEGGREST_1.44.1           statmod_1.5.0            
[127] igraph_2.0.3              memoise_2.0.1            
[129] DEoptimR_1.1-3            bit_4.5.0
```