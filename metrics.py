#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from sklearn import metrics


# In[ ]:

def clustering_evaluate(labels_true,labels_pred):
    
    NMI = metrics.normalized_mutual_info_score(labels_true, labels_pred)
    ARI = metrics.adjusted_rand_score(labels_true, labels_pred)
    print('NMI :',NMI)
    print('ARI :',ARI)
    return NMI,ARI
    
    
