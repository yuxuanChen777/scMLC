# scMLC
scMLC: an accurate and robust multiplex community detection method for single-cell multi-omics data
Yuxuan Chen, Ruiqing Zheng, Jin Liu, Min Li * 2023-03-22

Introduction
Clustering cells based on single-cell multi-modal sequencing technologies provides an unprecedented opportunity to create the high-resolution cell atlas, reveal cellular critical states and study health and disease. However, it is still a challenging task to effectively integrate different sequencing data for cell clustering. Motivated by the successful application of Louvain in scRNA-seq data, we propose a single-cell multi-modal Louvain clustering framework, called scMLC, to tackle the problem.

Workflow


Overview
The first step is feature selection;

The second step is construction of single-modal and cross-modal cell-to-cell networks;

The third step is weighting cells;

The last step is clustering by the multiplex network Louvain.

Usage
import imputation.py into bayesmodel.py
run the bayesmodel.py
