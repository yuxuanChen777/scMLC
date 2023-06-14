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

# Usage

Run the code in this way
```
run the lymph.R  
run the Louvain.ipynb
```
preprocessing
```
getwd()
setwd('/home/yxchen/R/R_data/Seurat/Lymph')
rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)
library(Signac)
library(RGCCA)
library(scry)
#
#


library(foreach)
library(doParallel)
library(RcppAnnoy)

getDoParWorkers( )    #查看注册了多少个核，配合doMC package中的registerDoMC( )使用
getDoParRegistered( ) # 查看doPar是否注册；如果没有注册返回FALSE
getDoParName( )       #查看已经注册的doPar的名字
getDoParVersion( )    #查看已经注册的doPar的version
getDoParWorkers()
detectCores() #
registerDoParallel(30)

mix_rna<-read.table('lymph_rna.tsv',sep='\t') #input data
mix_atac<-read.table('lymph_atac.tsv',sep='\t')


# #############################################################################################
rna <- CreateSeuratObject(counts = mix_rna, project = "mix_rna",min.cells = 10)
atac <- CreateSeuratObject(counts=mix_atac,project = 'mix_atac',assay='ATAC',min.cells = 10)
############################################################################################
#                       第二层scry特征提取的RNA网络                                        #                         
############################################################################################
rna <- SCTransform(rna,method = "glmGamPoi",verbose = FALSE) %>% RunPCA()
ElbowPlot(rna)
rna <- FindNeighbors(rna, dims = 1:10)#25
rna_layer <- as.data.frame(rna@graphs$SCT_nn)
write.table(rna_layer,file = '/home/yxchen/JupyterNotebook/Python_Louvain_op/PY_lymph_rna_layer.csv',sep = ',',row.names = T,col.names = T) #Output scRNA-seq single-mode network
############################################################################################
#                       第三层scry特征提取的ATAC网络                                       #                         
############################################################################################
atac_tmp <-devianceFeatureSelection(as.matrix(mix_atac))
atac_tmp <- as.data.frame(atac_tmp)
atac_tmp$peak_name <- rownames(atac_tmp)
atac_tmp <- atac_tmp[order(-atac_tmp$atac_tmp),]
atac_var <- atac_tmp[1:35000,]

atac_var <- atac_var$peak_name
atac@assays$ATAC@var.features <- atac_var
atac <- RunTFIDF(atac)
atac <- RunSVD(atac)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 2:10)
atac_layer <- as.data.frame(atac@graphs$ATAC_nn)
write.table(atac_layer,file = '/home/yxchen/JupyterNotebook/Python_Louvain_op/PY_lymph_atac_layer.csv',sep = ',',row.names = T,col.names = T)#Output scATAC-seq single-mode network

############################################################################################
#         第一层CCA特征提取网络（此处获得的是RNA，ATAC分别的两层，稍后会网络融合）         #
############################################################################################
rna <- CreateSeuratObject(counts = mix_rna, project = "mix_rna",min.cells = 10)
atac <- CreateSeuratObject(counts=mix_atac,project = 'mix_atac',assay='ATAC',min.cells = 10)

rna <- SCTransform(rna,method = "glmGamPoi",verbose = FALSE)
###########################################################################################
atac_tmp <-devianceFeatureSelection(as.matrix(mix_atac))
atac_tmp <- as.data.frame(atac_tmp)
atac_tmp$peak_name <- rownames(atac_tmp)
atac_tmp <- atac_tmp[order(-atac_tmp$atac_tmp),]
atac_var <- atac_tmp[1:30000,]
atac_var <- atac_var$peak_name
atac@assays$ATAC@var.features <- atac_var
atac <- RunTFIDF(atac)
atac <- ScaleData(atac,do.center = TRUE,do.scale = TRUE)

#########################################################################################
rna_m <-  as.matrix(rna@assays$SCT@scale.data)
atac_m <- as.matrix(atac@assays$ATAC@scale.data)

rna_f <- t(rna_m)
atac_f <-t(atac_m)
A <- list(rna_f,atac_f)
# sapply(A, head)
C = matrix(c(0, 1, 1, 0), 2, 2) 

start<-Sys.time()
sgcca.sci = sgcca(A, C,c1=c(1,0.5),ncomp = c(1, 1),scheme = "centroid",verbose = FALSE)# SGCCA:Cross-modal feature selection
end<-Sys.time()
runningtime<-end-start
print(runningtime) #1.43699 hours

nb_zero_rna = sum(sgcca.sci$a[[1]] != 0)
nb_zero_atac = sum(sgcca.sci$a[[2]] != 0)

rna_ff <- as.data.frame(sgcca.sci$a[[1]])
rna_ff$V2 <- rownames(rna_ff)
rna_ff <- rna_ff[rna_ff$V1!=0,]
rna_cca <- rna_ff$V2

atac_ff <- as.data.frame(sgcca.sci$a[[2]])
atac_ff$V2 <- rownames(atac_ff)
atac_ff <- atac_ff[atac_ff$V1!=0,]
atac_cca <- atac_ff$V2
#######################################################################
atac@assays$ATAC@var.features <- atac_cca
rna@assays$RNA@var.features <-rna_cca
#######################################################################
# #######################################################################
# #                      构建新的点权方法                               #
# #######################################################################
rna_mm <- as.data.frame(rna@assays$SCT@scale.data)
atac_mm <-as.data.frame(atac@assays$ATAC@scale.data)

K=20
F_num=nrow(rna_mm)
RNA_annoy <- new(AnnoyEuclidean,F_num)#计算的是欧氏距离
for (i in 1:ncol(rna_mm))
{
  RNA_annoy$addItem(i,rna_mm[,i])
}

RNA_annoy$build(50)

F_num_a=nrow(atac_mm)
ATAC_annoy <- new(AnnoyEuclidean,F_num_a)#计算的是欧氏距离
for (i in 1:ncol(atac_mm))
{
  ATAC_annoy$addItem(i,atac_mm[,i])
}

ATAC_annoy$build(50)


cell_weight <- data.frame(colnames(rna_mm))


res <- foreach (i = 1:ncol(rna_mm) ,.inorder = TRUE,.combine = 'c') %dopar% 
  {
    res_list <- RNA_annoy$getNNsByItemList(i,K,-1,TRUE)
    res_list$distance <- res_list$distance/res_list$distance[5] 
    A <- res_list$item
    mean_rna <- sum(res_list$distance)/20
    
    
    
    res_list_atac <- ATAC_annoy$getNNsByItemList(i,K,-1,TRUE)
    res_list_atac$distance <- res_list_atac$distance/res_list_atac$distance[5]
    B <- res_list_atac$item
    mean_atac <- sum(res_list_atac$distance)/20
    
    
    C <- intersect(A,B)
    D <- union(A,B)
    w <- (length(D)/length(C))
    # w <- (mean_rna+mean_atac)*(length(D)/length(C))
    return (1/w)
  }

cell_weight$W <-res


write.table(cell_weight,file='/home/yxchen/JupyterNotebook/Python_Louvain_op/lymph_cell_weight.csv',sep = ',',col.names = T)
#Output the weights of each node in the guidance layer network
#################################################################################################################################
rna <- RunPCA(rna,features = rna_cca)
atac <- RunSVD(atac,features = atac_cca)
rna <- FindNeighbors(rna,dims=1:10)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 2:10)
rna_knn <- as.data.frame(rna@graphs$SCT_nn)
atac_knn <- as.data.frame(atac@graphs$ATAC_nn)

write.table(rna_knn,file = '/home/yxchen/JupyterNotebook/Python_Louvain_op/PY_lymph_RNA_SGCCA_layer.csv',sep = ',',row.names = T,col.names = T) #Output scRNA monolayer network after cross-modal feature selection
write.table(atac_knn,file = '/home/yxchen/JupyterNotebook/Python_Louvain_op/PY_lymph_ATAC_SGCCA_layer.csv',sep = ',',row.names = T,col.names = T) #Output scATAC monolayer network after cross-modal feature selection

```
