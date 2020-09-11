#!/usr/bin/env python
# coding: utf-8

# In[59]:


# import libraries and packages necessary for analysis

# Load the Libraries

# import some packages that may be helpful for us

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


# In[60]:


# define a random seed to replicate the results
RSEED = 50


# In[72]:


# Read in the data set




X= pd.read_csv("~/storage/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/TCGA-BRCA_metastatic_data_RNAseq.csv")
X.drop(X.columns[0], axis=1)
# Look at the shape of the data to construct the input layer. 
print('We have {} instances of data with {} variables'.format(*X.shape))


# In[73]:


y = X[["Bone","Lung","Liver"]].to_numpy()
n_classes = y.shape[1]
patients = np.array(X.pop('barcode'))
# remove the features labels
X = X.drop(["Bone","Lung","Liver"], axis = 1) 


# In[74]:


from sklearn.metrics import accuracy_score
import numpy as np
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score


# In[75]:


# Labels are the values we want to predict, retrun and pop them off after
# List of features for later use
feature_list = list(X.columns)
# Convert to numpy array
features = np.array(X)


# In[86]:


train, test, train_labels, test_labels = train_test_split(X,
                                         y,
                                         stratify =y,
                                         test_size = 0.3, 
                                         random_state = RSEED)

model = RandomForestClassifier(n_estimators=1000, 
                               random_state=RSEED, 
                               max_features = 'sqrt',
                               n_jobs=-1, verbose = 1)


# In[87]:


model.fit(train, train_labels)
# Training predictions (to demonstrate overfitting)
train_rf_predictions = model.predict(train)
train_rf_probs = model.predict_proba(train)

# Testing predictions (to determine performance)
rf_predictions = model.predict(test)
rf_probs = model.predict_proba(test)


# In[88]:


model.score(train,train_labels)
#model.predict_proba(test)


# In[89]:


y_score =model.score(test,test_labels)


# In[90]:


from sklearn.metrics import confusion_matrix
import itertools
cm = confusion_matrix(
    test_labels.argmax(axis=1), rf_predictions.argmax(axis=1))


# In[81]:


y_score


# In[ ]:


# test accuracy for each class

