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

Please refer to the vignettes. All data required for vignettes has been uploaded to the onedrive, which can be access by following link: https://1drv.ms/f/c/60d1e294ce61ccdc/EoF_83k6p9tAvJYrm5BSIgYBd27wGIrUtfdNL4s-5mk8EA?e=ugjUeZ

### Required Package

In brief, you need to install the following packages: 

1. R packages:

    ENIGMA, abind, data.table, pcalg, reticulate, scater, WGCNA, coop, reticulate, stats, xgboost, zellkonverter
    
2. Python packages:

    torch, numpy, pandas, matplotlib, seaborn, anndata, sklearn, qnorm, combat

We strongly recommend you to install the pcalg we provided in the Github repository, as we have revise some of the functions provided in the orignal package. Details can be referred to our previously published nature machine intelligence paper.

Additionally, we attached the output of sessionInfo() to show the exact version of the package we used in our analysis.

```
> sessionInfo()
R version 4.5.1 (2025-06-13)
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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] reticulate_1.44.0           CausalCellInfer_1.0        
 [3] scater_1.38.0               ggplot2_4.0.0              
 [5] scuttle_1.20.0              SingleCellExperiment_1.32.0
 [7] SummarizedExperiment_1.40.0 Biobase_2.70.0             
 [9] GenomicRanges_1.62.0        Seqinfo_1.0.0              
[11] IRanges_2.44.0              S4Vectors_0.48.0           
[13] BiocGenerics_0.56.0         generics_0.1.4             
[15] MatrixGenerics_1.22.0       matrixStats_1.5.0          
[17] zellkonverter_1.20.0        stringr_1.6.0              
[19] data.table_1.17.8          

loaded via a namespace (and not attached):
  [1] splines_4.5.1                      filelock_1.0.3                    
  [3] tibble_3.3.0                       preprocessCore_1.72.0             
  [5] graph_1.88.0                       XML_3.99-0.20                     
  [7] rpart_4.1.24                       lifecycle_1.0.4                   
  [9] CVST_0.2-3                         fastcluster_1.3.0                 
 [11] edgeR_4.8.0                        doParallel_1.0.17                 
 [13] lattice_0.22-7                     MASS_7.3-65                       
 [15] backports_1.5.0                    magrittr_2.0.4                    
 [17] limma_3.66.0                       Hmisc_5.2-4                       
 [19] rmarkdown_2.30                     DBI_1.2.3                         
 [21] RColorBrewer_1.1-3                 abind_1.4-8                       
 [23] sfsmisc_1.1-22                     purrr_1.2.0                       
 [25] nnet_7.3-20                        sva_3.58.0                        
 [27] ggrepel_0.9.6                      irlba_2.3.5.1                     
 [29] ENIGMA_0.1.5                       genefilter_1.92.0                 
 [31] annotate_1.88.0                    codetools_0.2-20                  
 [33] DelayedArray_0.36.0                tidyselect_1.2.1                  
 [35] UCSC.utils_1.6.0                   farver_2.1.2                      
 [37] ScaledMatrix_1.18.0                viridis_0.6.5                     
 [39] dynamicTreeCut_1.63-1              base64enc_0.1-3                   
 [41] jsonlite_2.0.0                     BiocNeighbors_2.4.0               
 [43] e1071_1.7-16                       Formula_1.2-5                     
 [45] survival_3.8-3                     iterators_1.0.14                  
 [47] foreach_1.5.2                      tools_4.5.1                       
 [49] Rcpp_1.1.0                         glue_1.8.0                        
 [51] gridExtra_2.3                      SparseArray_1.10.1                
 [53] xfun_0.54                          mgcv_1.9-4                        
 [55] GenomeInfoDb_1.46.0                dplyr_1.1.4                       
 [57] withr_3.0.2                        BiocManager_1.30.26               
 [59] fastmap_1.2.0                      basilisk_1.22.0                   
 [61] GeneralisedCovarianceMeasure_0.2.0 digest_0.6.37                     
 [63] rsvd_1.0.5                         R6_2.6.1                          
 [65] colorspace_2.1-2                   GO.db_3.22.0                      
 [67] RSQLite_2.4.3                      tidyr_1.3.1                       
 [69] pcalg_2.7-12                       corpcor_1.6.10                    
 [71] robustbase_0.99-6                  class_7.3-23                      
 [73] httr_1.4.7                         htmlwidgets_1.6.4                 
 [75] S4Arrays_1.10.0                    pkgconfig_2.0.3                   
 [77] gtable_0.3.6                       blob_1.2.4                        
 [79] S7_0.2.0                           impute_1.84.0                     
 [81] XVector_0.50.0                     htmltools_0.5.8.1                 
 [83] RBGL_1.86.0                        clue_0.3-66                       
 [85] scales_1.4.0                       png_0.1-8                         
 [87] knitr_1.50                         rstudioapi_0.17.1                 
 [89] ggm_2.5.2                          reshape2_1.4.4                    
 [91] checkmate_2.3.3                    nlme_3.1-168                      
 [93] bdsmatrix_1.3-7                    proxy_0.4-27                      
 [95] cachem_1.1.0                       parallel_4.5.1                    
 [97] vipor_0.4.7                        foreign_0.8-90                    
 [99] AnnotationDbi_1.72.0               pillar_1.11.1                     
[101] grid_4.5.1                         fastICA_1.2-7                     
[103] vctrs_0.6.5                        BiocSingular_1.26.0               
[105] beachmat_2.26.0                    xtable_1.8-4                      
[107] cluster_2.1.8.1                    beeswarm_0.4.0                    
[109] htmlTable_2.4.3                    evaluate_1.0.5                    
[111] cli_3.6.5                          locfit_1.5-9.12                   
[113] compiler_4.5.1                     rlang_1.1.6                       
[115] crayon_1.5.3                       plyr_1.8.9                        
[117] ggbeeswarm_0.7.2                   stringi_1.8.7                     
[119] viridisLite_0.4.2                  WGCNA_1.73                        
[121] BiocParallel_1.44.0                nnls_1.6                          
[123] Biostrings_2.78.0                  coop_0.6-3                        
[125] Matrix_1.7-4                       dir.expiry_1.18.0                 
[127] bit64_4.6.0-1                      KEGGREST_1.50.0                   
[129] statmod_1.5.1                      kernlab_0.9-33                    
[131] igraph_2.2.1                       memoise_2.0.1                     
[133] DEoptimR_1.1-4                     bit_4.6.0                         
[135] xgboost_1.7.11.1           
```