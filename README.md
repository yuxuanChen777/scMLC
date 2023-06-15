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
**Preprocessing lymph.R**
```
mix_rna<-read.table('lymph_rna.tsv',sep='\t') #input data
mix_atac<-read.table('lymph_atac.tsv',sep='\t')
```
**multiplex network community detection  Louvain.ipynb**
```
cell_weight_df = pd.read_csv('lymph_cell_weight.csv',sep=',',header=0) # cell_weight lymph
second_layer=pd.read_csv('PY_lymph_rna_layer.csv' ,header=0,index_col=0,sep=',')# scRNA-seq single-mode network  
df01_1=pd.read_csv('PY_lymph_RNA_SGCCA_layer.csv',header=0,index_col=0,sep=',') # scRNA-seq cross-mode network  
df01_2=pd.read_csv('PY_lymph_ATAC_SGCCA_layer.csv',header=0,index_col=0,sep=',') # scATAC-seq cross-mode network  
third_layer=pd.read_csv('PY_lymph_atac_layer.csv',header=0,index_col=0,sep=',') # scATAC-seq single-mode network  
RES.to_csv('PY_lymph_Tri_no_w.tsv',sep='\t',header=None,index=False) # multi-plex network  
nodes, layer_data = load_multilayer_network('PY_lymph_Tri_no_w.tsv') # load multi-plex network
label=pd.read_csv('lymph_cluster.csv',sep = ',',header=0) #True Label
```

