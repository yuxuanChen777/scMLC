#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from sklearn import metrics


# In[ ]:

# def accuracy(labels_true, labels_pred):
#     clusters = np.unique(labels_pred)
#     labels_true = np.reshape(labels_true, (-1, 1))
#     labels_pred = np.reshape(labels_pred, (-1, 1))
#     count = []
#     for c in clusters:
#         idx = np.where(labels_pred == c)[0]
#         labels_tmp = labels_true[idx, :].reshape(-1)
#         count.append(np.bincount(labels_tmp).max())
#     return np.sum(count) / labels_true.shape[0]

# def get_rand_index_and_f_measure(labels_true, labels_pred, beta=1.):
#     (tn, fp), (fn, tp) = pair_confusion_matrix(labels_true, labels_pred)
#     ri = (tp + tn) / (tp + tn + fp + fn)
#     ari = 2. * (tp * tn - fn * fp) / ((tp + fn) * (fn + tn) + (tp + fp) * (fp + tn))
#     p, r = tp / (tp + fp), tp / (tp + fn)
#     f_beta = (1 + beta**2) * (p * r / ((beta ** 2) * p + r))
#     return f_beta

def clustering_evaluate(labels_true,labels_pred):
    
    NMI = metrics.normalized_mutual_info_score(labels_true, labels_pred)
    ARI = metrics.adjusted_rand_score(labels_true, labels_pred)
    FM = metrics.fowlkes_mallows_score(labels_true, labels_pred)
    print('NMI :',NMI)
    #print('v_mea :',v_mea)
    print('ARI :',ARI)
    #print('FM :',FM)
    #print('Purity :',Purity)
    return NMI,ARI
    
    
