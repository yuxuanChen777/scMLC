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

**Run the code in this way**
```
run the lymph.R  
run the Louvain.ipynb
```
**Preprocessing**
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

library(foreach)
library(doParallel)
library(RcppAnnoy)

getDoParWorkers( )  
getDoParRegistered( ) 
getDoParName( )       
getDoParVersion( )   
getDoParWorkers()
detectCores() 
registerDoParallel(30)

mix_rna<-read.table('lymph_rna.tsv',sep='\t') #input data
mix_atac<-read.table('lymph_atac.tsv',sep='\t')


# ##########################################################################################
rna <- CreateSeuratObject(counts = mix_rna, project = "mix_rna",min.cells = 10)
atac <- CreateSeuratObject(counts=mix_atac,project = 'mix_atac',assay='ATAC',min.cells = 10)
############################################################################################
#                      scRNA-seq single-mode network                                       #                         
############################################################################################
rna <- SCTransform(rna,method = "glmGamPoi",verbose = FALSE) %>% RunPCA()
ElbowPlot(rna)
rna <- FindNeighbors(rna, dims = 1:10)#25
rna_layer <- as.data.frame(rna@graphs$SCT_nn)
write.table(rna_layer,file = '/home/yxchen/JupyterNotebook/Python_Louvain_op/PY_lymph_rna_layer.csv',sep = ',',row.names = T,col.names = T) #Output scRNA-seq single-mode network
############################################################################################
#                       scATAC-seq single-mode network                                     #                         
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
#                              Cross-modal feature selection                               #
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
# #     Assigns weights to nodes in the guidance layer network          #
# #######################################################################
rna_mm <- as.data.frame(rna@assays$SCT@scale.data)
atac_mm <-as.data.frame(atac@assays$ATAC@scale.data)

K=20
F_num=nrow(rna_mm)
RNA_annoy <- new(AnnoyEuclidean,F_num)
for (i in 1:ncol(rna_mm))
{
  RNA_annoy$addItem(i,rna_mm[,i])
}

RNA_annoy$build(50)

F_num_a=nrow(atac_mm)
ATAC_annoy <- new(AnnoyEuclidean,F_num_a)
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
###multiplex network community detection###
```
import networkx as nx
import random as rd
import math
import pandas as pd
from community_louvain import *
import operator
import metrics
import umap
import numpy as np
# cell_weight_df = pd.read_csv('sci_CAR_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('P0_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('LLL_CTRL_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('Paired_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('celllinemixture_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('PBMC10_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('NEAT_seq_MIX_cell_weight.csv',sep=',',header=0)
cell_weight_df = pd.read_csv('lymph_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('NEAT_seq_CD4T_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('Ad_cell_weight.csv',sep=',',header=0)
# cell_weight_df = pd.read_csv('PBMC3_cell_weight.csv',sep=',',header=0)
cell_weight_df.columns = ['Cell','Weight']
# cell_weight_df.sort_values("Weight",inplace=True, ascending=True)
cell_weight_df.sort_values("Weight",inplace=True, ascending=False)
Cell_sort = cell_weight_df['Cell'].tolist()
len(Cell_sort)
7058
def load_multilayer_network(input_file):
    layer_data = {}#key为0，1对应的value是边
    nodes_lst = []#存的是全部的node
    nodes = {}#字典 将细胞的名字对应为数字
    cnt=0
    with open(input_file, 'r') as _file:
        for line in _file:
            cnt = cnt + 1
            line = line.strip()
            data = [x.strip() for x in line.split('\t')]
            nodes_lst.append(data[0])
            nodes_lst.append(data[1])
    _file.close()

#     print('type(nodes_lst) :',type(nodes_lst))
    nodes_lst = list(set(nodes_lst))
#     print('type(nodes_lst) :',type(nodes_lst))
    for i in range(len(nodes_lst)):
        nodes[nodes_lst[i]] = i
    
    with open(input_file, 'r') as _file:
        for line in _file:
            line = line.strip()
            data = [x.strip() for x in line.split('\t')]
            if int(float(data[3])) not in layer_data.keys():
                layer_data[int(float(data[3]))] = []
    
            layer_data[int(float(data[3]))].append((nodes[data[0]], nodes[data[1]],float(data[2])))        
    _file.close()
    
#     print(nodes)
#     print(layer_data)
    return nodes, layer_data
def make_multilayer(layer_data, node_lst):
    multilayer_network = {} #字典
#     import pdb
#     pdb.set_trace()
    for i in layer_data.keys():
            multilayer_network[i] = nx.Graph()
            for edge in layer_data[i]:
                #print('edge :',edge)
                multilayer_network[i].add_weighted_edges_from([(edge[0], edge[1],edge[2])])
                #print('type(edge[2]):',type(edge[2]))
            multilayer_network[i].add_nodes_from(node_lst)
    return multilayer_network


# second_layer=pd.read_csv('PY_sci_CAR_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_An_Ultra_Paired_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_cellmixture_SNARE_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_P0_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_NEAT_seq_MIX_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_PBMC10_scMVP_rna_layer.csv' ,header=0,index_col=0,sep=',')
second_layer=pd.read_csv('PY_lymph_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_Ad_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_PBMC3_rna_layer.csv' ,header=0,index_col=0,sep=',')
# second_layer=pd.read_csv('PY_LLL_CTRL_rna_layer.csv' ,header=0,index_col=0,sep=',')

# second_layer=pd.read_csv('PY_NEAT_seq_CD4T_rna_layer.csv' ,header=0,index_col=0,sep=',')
df_triple=second_layer.stack().reset_index()
df_triple.columns = ['U','V','Value']
df_triple_clear=df_triple[df_triple.Value!=0]
df_triple_clear_Final=df_triple_clear[df_triple_clear['U']!=df_triple_clear['V']]
df_triple_clear_Final['Layer']= [1 for index in range(len(df_triple_clear_Final))]

# df01_1=pd.read_csv('PY_Ad_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_Ad_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_An_Ultra_Paired_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_An_Ultra_Paired_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
df01_1=pd.read_csv('PY_lymph_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
df01_2=pd.read_csv('PY_lymph_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_sci_CAR_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_sci_CAR_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_cellmixture_SNARE_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_cellmixture_SNARE_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_P0_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_P0_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_PBMC10_scMVP_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_PBMC10_scMVP_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_NEAT_seq_MIX_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_NEAT_seq_MIX_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_3=pd.read_csv('PY_NEAT_seq_MIX_ADT_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_NEAT_seq_CD4T_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_NEAT_seq_CD4T_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_3=pd.read_csv('PY_NEAT_seq_CD4T_ADT_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_PBMC3_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_PBMC3_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_1=pd.read_csv('PY_LLL_CTRL_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',')
# df01_2=pd.read_csv('PY_LLL_CTRL_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',')
#

first_layer=df01_1+df01_2
# first_layer=df01_1+df01_2+df01_3
df01_triple=first_layer.stack().reset_index()
df01_triple.columns = ['U','V','Value']
df01_triple_clear=df01_triple[df01_triple.Value!=0]
df01_triple_clear_Final=df01_triple_clear[df01_triple_clear['U']!=df01_triple_clear['V']]
df01_triple_clear_Final['Layer']= [0 for index in range(len(df01_triple_clear_Final))]
df01_triple_clear_Final['Value']=df01_triple_clear_Final['Value'].map(lambda x : 1)

# third_layer=pd.read_csv('PY_sci_CAR_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_An_Ultra_Paired_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_cellmixture_SNARE_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_P0_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_PBMC10_scMVP_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_NEAT_seq_MIX_atac_layer.csv',header=0,index_col=0,sep=',')
third_layer=pd.read_csv('PY_lymph_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_Ad_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_PBMC3_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_NEAT_seq_CD4T_atac_layer.csv',header=0,index_col=0,sep=',')
# third_layer=pd.read_csv('PY_LLL_CTRL_atac_layer.csv',header=0,index_col=0,sep=',')
df02_triple=third_layer.stack().reset_index()
df02_triple.columns = ['U','V','Value']
df02_triple_clear=df02_triple[df02_triple.Value!=0]
df02_triple_clear_Final=df02_triple_clear[df02_triple_clear['U']!=df02_triple_clear['V']]
df02_triple_clear_Final['Layer']= [2 for index in range(len(df02_triple_clear_Final))]
df02_triple_clear_Final

# fourth_layer=pd.read_csv('PY_LLL_CTRL_adt_layer.csv',header=0,index_col=0,sep=',')
# fourth_layer=pd.read_csv('PY_NEAT_seq_MIX_adt_layer.csv',header=0,index_col=0,sep=',')
# # # fourth_layer=pd.read_csv('PY_NEAT_seq_CD4T_adt_layer.csv',header=0,index_col=0,sep=',')
# df03_triple=fourth_layer.stack().reset_index()
# df03_triple.columns = ['U','V','Value']
# df03_triple_clear=df03_triple[df03_triple.Value!=0]
# df03_triple_clear_Final=df03_triple_clear[df03_triple_clear['U']!=df03_triple_clear['V']]
# df03_triple_clear_Final['Layer']= [3 for index in range(len(df03_triple_clear_Final))]
# df03_triple_clear_Final

# # frames=[df01_triple_clear_Final]
# frames=[df01_triple_clear_Final,df_triple_clear_Final,df02_triple_clear_Final,df03_triple_clear_Final]
# RES=pd.concat(frames)
# RES.to_csv('PY_NEAT_seq_MIX_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_NEAT_seq_MIX_Tri_no_w.tsv',sep='\t',header=None,index=False)
# # RES.to_csv('PY_NEAT_seq_CD4T_Tri_no_w.tsv',sep='\t',header=None,index=False)

frames=[df01_triple_clear_Final,df_triple_clear_Final,df02_triple_clear_Final]
RES=pd.concat(frames)
# # RES.to_csv('PY_Ad_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_NEAT_seq_MIX_Tri_no_w.tsv',sep='\t',header=None,index=False)
# # RES.to_csv('PY_NEAT_seq_CD4T_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_PBMC10_scMVP_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_sci_CAR_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_An_Ultra_Paired_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_cellmixture_SNARE_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_P0_Tri_no_w.tsv',sep='\t',header=None,index=False)
RES.to_csv('PY_lymph_Tri_no_w.tsv',sep='\t',header=None,index=False)
# # RES.to_csv('PY_PBMC3_Tri_no_w.tsv',sep='\t',header=None,index=False)
# RES.to_csv('PY_LLL_CTRL_Four_no_w.tsv',sep='\t',header=None,index=False)
# nodes, layer_data = load_multilayer_network('PY_NEAT_seq_MIX_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_NEAT_seq_CD4T_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_cellmixture_SNARE_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_P0_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_LLL_CTRL_Four_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_An_Ultra_Paired_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_PBMC10_scMVP_Tri_no_w.tsv')
nodes, layer_data = load_multilayer_network('PY_lymph_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_An_Ultra_Paired_rna_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_Ad_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_PBMC3_Tri_no_w.tsv')
# nodes, layer_data = load_multilayer_network('PY_sci_CAR_Tri_no_w.tsv')
nodes_lst = sorted(list(nodes.values()))
len(layer_data)
4
len(nodes)
7058
Cell_Sort_Id=list()
for i in range(len(Cell_sort)):
    Cell_Sort_Id.append(nodes[Cell_sort[i]])
cell_weight_df['ID']=Cell_Sort_Id
Cell_Id_W_dict = dict(zip(cell_weight_df['ID'].tolist(),cell_weight_df['Weight'].tolist()))
len(Cell_Id_W_dict)


multilayer_network=make_multilayer(layer_data, nodes_lst) #multilayer_network是字典，对应的values是值
7058
partition = best_partition(multilayer_network,Cell_Id_W_dict,weight='weight',resolution=0.9)#3 0.433(P0)0.9 0.87(lymph) 0.508(PBMC10)
community_num = len(set(partition.values()))
print('Community_Num :',community_num)
resolution : 2.85
_one_level_first
_one_level_second
_one_level_second
Community_Num : 6
names = {}
for n in nodes:
#     print('n :',n)
#     print('nodes[n] :',nodes[n])
    names[nodes[n]] = n
B=[]
C=[]

for j in partition.keys():
        B.append(names[j])
for k in partition.values():
        C.append(k)
Res=dict(zip(B, C))
Res = dict(sorted(Res.items(), key=operator.itemgetter(0)))  #按照key值升序
# print('Res :',Res)
# Res
# jsObj = json.dumps(Res)
 
# # fileObject = open('lymph_C_Random.json', 'w')
# fileObject = open('lymph_C.json', 'w')
# fileObject.write(jsObj)
# fileObject.close()
# label=pd.read_csv('celllinemixture_cluster_omit.csv',sep = ',',header=0)
# label=pd.read_csv('An_Ultra_cluster_omit_Final_intersect_dc.csv',sep=',',header=0)
# label=pd.read_csv('sci_CAR_cell_atac.tsv',sep='\t',header=0)
# label=pd.read_csv('HT_Ad_cluster.csv',sep = ',',header=0)
# label=pd.read_csv('HT_P0_cluster.csv',sep = ',',header=0) 
# label=pd.read_csv('PBMC10_scMVP_cluster.csv',sep = ',',header=0)
# label=pd.read_csv('NEAT_seq_MIX_cluster_new.csv',sep = ',',header=0)
label=pd.read_csv('lymph_cluster.csv',sep = ',',header=0)
# label=pd.read_csv('NEAT_seq_CD4_T_cluster.csv',sep = ',',header=0)
# label=pd.read_csv('SNARE_cellline_cluster.csv',sep = ',',header=0)
# label=pd.read_csv('PBMC3_cluster.csv',sep = ',',header=0)
# label=pd.read_csv('NEAT_seq_MIX_cluster_Three.csv',sep = ',',header=0)
# label=pd.read_csv('LLL_CTRL_cluster.csv',sep = ',',header=0)
# label_index=label.sort_values(by="cell_id",inplace=False)
# label_index
# Y=np.array(label_index['cluster_id'])  
# len(set(Y))
# Y=np.array(label_index['cluster_id'])  
# len(set(Y))
# #An_Ultra
# Y=np.array(label['Ident_ID'])
# len(set(Y))
# # # # # label_index=label.sort_index()  #HT_celllinemixture
# Y=np.array(label['A_Type'])
# # #ATAC_Y=np.array(label['A_Type'])
# print('len(Y) :',len(set(Y)))
# Y=np.array(label['Ident'])  #P0 Ad
# len(set(Y))
Y=np.array(label['cluster_id'])  #PBMC10_scMVP  lymph SHARE sci_CAR NEAT_seq_MIX
len(set(Y))
set(Y)
{1, 2, 3, 4, 5, 6}
Pre=[]
for k in Res.values():
        Pre.append(k)
Pre=np.array(Pre)
print(len(set(Pre)))
metrics.clustering_evaluate(Y,Pre)
```

