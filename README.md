# scMLC
scMLC: an accurate and robust multiplex community detection method for single-cell multi-omics data

# Introduction
Clustering cells based on single-cell multi-modal sequencing technologies provides an unprecedented opportunity to create the high-resolution cell atlas, reveal cellular critical states and study health and disease. However, it is still a challenging task to effectively integrate different sequencing data for cell clustering. Motivated by the successful application of Louvain in scRNA-seq data, we propose a single-cell multi-modal Louvain clustering framework, called scMLC, to tackle the problem.  

# Datasets
Lymph: RNA + ATAC

# Workflow


# Overview
 * The first step is feature selection;  

 * The second step is construction of single-modal and cross-modal cell-to-cell networks;  

 * The third step is weighting cells;  

 * The last step is clustering by the multiplex network Louvain.

# R Packages version
```
sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RcppAnnoy_0.0.20            scry_1.8.0                  Signac_1.10.0               reticulate_1.30             patchwork_1.1.2            
 [6] SeuratObject_4.1.3          Seurat_4.3.0.1              RGCCA_2.1.2                 ACAT_0.91                   umap_0.2.10.0              
[11] ggh4x_0.2.4                 SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
[16] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
[21] matrixStats_1.0.0           scDesign3_0.99.5            lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0              
[26] dplyr_1.1.2                 purrr_1.0.1                 readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1               
[31] ggplot2_3.4.2               tidyverse_2.0.0             doParallel_1.0.17           iterators_1.0.14            foreach_1.5.2              
[36] MiXcan_0.1.0                Matrix_1.5-4.1             

loaded via a namespace (and not attached):
  [1] utf8_1.2.3                spatstat.explore_3.2-1    tidyselect_1.2.0          htmlwidgets_1.6.2         grid_4.2.0               
  [6] BiocParallel_1.32.6       Rtsne_0.16                munsell_0.5.0             ScaledMatrix_1.4.1        codetools_0.2-18         
 [11] ica_1.0-3                 future_1.32.0             miniUI_0.1.1.1            withr_2.5.0               spatstat.random_3.1-5    
 [16] colorspace_2.1-0          progressr_0.13.0          highr_0.10                knitr_1.43                rstudioapi_0.14          
 [21] ROCR_1.0-11               tensor_1.5                listenv_0.9.0             labeling_0.4.2            GenomeInfoDbData_1.2.9   
 [26] polyclip_1.10-4           farver_2.1.1              parallelly_1.36.0         vctrs_0.6.3               generics_0.1.3           
 [31] xfun_0.39                 timechange_0.2.0          R6_2.5.1                  rsvd_1.0.5                bitops_1.0-7             
 [36] spatstat.utils_3.0-3      DelayedArray_0.24.0       promises_1.2.0.1          scales_1.2.1              gtable_0.3.3             
 [41] beachmat_2.12.0           globals_0.16.2            goftest_1.2-3             rlang_1.1.1               RcppRoll_0.3.0           
 [46] splines_4.2.0             lazyeval_0.2.2            spatstat.geom_3.2-1       reshape2_1.4.4            abind_1.4-5              
 [51] httpuv_1.6.11             tools_4.2.0               ellipsis_0.3.2            gamlss.data_6.0-2         RColorBrewer_1.1-3       
 [56] ggridges_0.5.4            gamlss_5.4-12             Rcpp_1.0.10               plyr_1.8.8                sparseMatrixStats_1.8.0  
 [61] zlibbioc_1.44.0           RCurl_1.98-1.12           openssl_2.0.6             deldir_1.0-9              pbapply_1.7-0            
 [66] cowplot_1.1.1             zoo_1.8-12                ggrepel_0.9.3             cluster_2.1.3             magrittr_2.0.3           
 [71] data.table_1.14.8         RSpectra_0.16-1           glmGamPoi_1.8.0           scattermore_1.2           lmtest_0.9-40            
 [76] RANN_2.6.1                fitdistrplus_1.1-11       hms_1.1.3                 mime_0.12                 evaluate_0.21            
 [81] xtable_1.8-4              mclust_6.0.0              gridExtra_2.3             compiler_4.2.0            KernSmooth_2.23-20       
 [86] crayon_1.5.2              htmltools_0.5.5           later_1.3.1               tzdb_0.4.0                MASS_7.3-57              
 [91] cli_3.6.1                 igraph_1.4.3              pkgconfig_2.0.3           sp_1.6-1                  plotly_4.10.2            
 [96] spatstat.sparse_3.0-1     gamlss.dist_6.0-5         XVector_0.38.0            digest_0.6.31             sctransform_0.3.5        
[101] spatstat.data_3.0-1       Biostrings_2.64.1         leiden_0.4.3              fastmatch_1.1-3           uwot_0.1.14              
[106] DelayedMatrixStats_1.18.2 shiny_1.7.4               Rsamtools_2.12.0          lifecycle_1.0.3           nlme_3.1-162             
[111] jsonlite_1.8.5            viridisLite_0.4.2         askpass_1.1               fansi_1.0.4               pillar_1.9.0             
[116] lattice_0.20-45           fastmap_1.1.1             httr_1.4.6                survival_3.5-5            glue_1.6.2               
[121] png_0.1-8                 stringi_1.7.12            BiocSingular_1.12.0       irlba_2.3.5.1             future.apply_1.11.0 
```
# Python version
```
Python 3.7.11
```

# Usage

**Run the code in this way**
```
run the lymph.R  
run the Louvain.ipynb
```
**Preprocessing lymph.R**
```
mix_rna<-read.table('lymph_rna.tsv',sep='\t') #input data
mix_atac<-read.table('lymph_atac.tsv',sep='\t')
```
**Multiplex network community detection  Louvain.ipynb**
```
cell_weight_df = pd.read_csv('lymph_cell_weight.csv',sep=',',header=0) # cell_weight lymph
second_layer=pd.read_csv('PY_lymph_rna_layer.csv' ,header=0,index_col=0,sep=',')# scRNA-seq single-mode network  
df01_1=pd.read_csv('PY_lymph_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',') # scRNA-seq cross-mode network  
df01_2=pd.read_csv('PY_lymph_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',') # scATAC-seq cross-mode network  
third_layer=pd.read_csv('PY_lymph_atac_layer.csv',header=0,index_col=0,sep=',') # scATAC-seq single-mode network  
RES.to_csv('PY_lymph_Tri_no_w.tsv',sep='\t',header=None,index=False) # multiplex network  
nodes, layer_data = load_multilayer_network('PY_lymph_Tri_no_w.tsv') # load multiplex network
label=pd.read_csv('lymph_cluster.csv',sep = ',',header=0) #True Label
```

