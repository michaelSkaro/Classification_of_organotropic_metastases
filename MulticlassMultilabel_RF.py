#!/usr/bin/env python
# coding: utf-8

# In[178]:


# import libraries and packages necessary for analysis

# Load the Libraries

# import some packages that may be helpful for us

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


# In[179]:


# define a random seed to replicate the results
RSEED = 50


# In[180]:


# Read in the data set




features= pd.read_csv("~/storage/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/TCGA-BLCA_metastatic_data_RNAseq.csv")

# Look at the shape of the data to construct the input layer. 
print('We have {} instances of data with {} variables'.format(*features.shape))


# In[181]:


type(features)


# In[182]:


y = features[["Lymph Node","Prostate","Bladder","Lung","Bone","Liver","Pelvis"]].to_numpy()


# In[176]:


y


# In[184]:


# remove the features labels
features = features.drop(["Lymph Node","Prostate","Bladder","Lung","Bone","Liver","Pelvis"], axis = 1) 


# In[185]:


#labels are not balanced, we need to resample to balance them
from sklearn.utils import resample
from sklearn.metrics import accuracy_score
import numpy as np
from sklearn.multioutput import MultiOutputClassifier


# In[186]:


# Labels are the values we want to predict, retrun and pop them off after
# grab the patients data ids 
patients = np.array(features.pop('barcode'))
# Saving feature names for later use
feature_list = list(features.columns)
# Convert to numpy array
features = np.array(features)


# In[208]:



model = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, 
                               random_state=RSEED, 
                               max_features = 'sqrt',
                               n_jobs=-1, verbose = 1)).fit(features, y)


# In[210]:


model.estimators_


# In[ ]:


print("done")


# In[213]:


model.score(features,y)


# In[215]:


model.predict(features)


# In[ ]:


model.predict_proba(features)

