#!/usr/bin/env python
# coding: utf-8

# In[269]:


# import libraries and packages necessary for analysis

# Load the Libraries

# import some packages that may be helpful for us

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


# In[270]:


# define a random seed to replicate the results
RSEED = 50


# In[271]:


# Read in the data set




features= pd.read_csv("~/storage/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/TCGA-BRCA_metastatic_data_RNAseq.csv")

# Look at the shape of the data to construct the input layer. 
print('We have {} instances of data with {} variables'.format(*features.shape))


# In[272]:


type(features)


# In[273]:


y = features[["Bone","Lung","Liver"]].to_numpy()


# In[274]:


y


# In[275]:


# remove the features labels
features = features.drop(["Bone","Lung","Liver"], axis = 1) 


# In[276]:


#labels are not balanced, we need to resample to balance them
from sklearn.utils import resample
from sklearn.metrics import accuracy_score
import numpy as np
from sklearn.multioutput import MultiOutputClassifier


# In[277]:


# Labels are the values we want to predict, retrun and pop them off after
# grab the patients data ids 
patients = np.array(features.pop('barcode'))
# Saving feature names for later use
feature_list = list(features.columns)
# Convert to numpy array
features = np.array(features)


# In[278]:


# 30% examples in test data
train, test, train_labels, test_labels = train_test_split(features,
                                         y,
                                         stratify =y,
                                         test_size = 0.3, 
                                         random_state = RSEED)


# In[279]:


model = MultiOutputClassifier(RandomForestClassifier(n_estimators=100, 
                               random_state=RSEED, 
                               max_features = 'sqrt',
                               n_jobs=-1, verbose = 1))


# In[280]:


model.fit(train, train_labels)


# In[281]:


# Training predictions (to demonstrate overfitting)
train_rf_predictions = model.predict(train)
train_rf_probs = model.predict_proba(train)

# Testing predictions (to determine performance)
rf_predictions = model.predict(test)
rf_probs = model.predict_proba(test)


# In[282]:


model.estimators_


# In[283]:


model.score(features,y)


# In[285]:


model.predict(features)


# In[286]:


model.predict_proba(features)


# In[287]:


rf_probs


# In[288]:


model.score(train,train_labels)


# In[289]:


model.score(test,test_labels)


# In[290]:


train_rf_predictions
train_rf_probs

# Testing predictions (to determine performance)
rf_predictions
rf_probs


# In[292]:


from sklearn.metrics import confusion_matrix
import itertools

def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Oranges):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    Source: http://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    # Plot the confusion matrix
    plt.figure(figsize = (10, 10))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title, size = 24)
    plt.colorbar(aspect=4)
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, size = 14)
    plt.yticks(tick_marks, classes, size = 14)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    
    # Labeling the plot
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt), fontsize = 20,
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")
        
    plt.grid(None)
    plt.tight_layout()
    plt.ylabel('True label', size = 18)
    plt.xlabel('Predicted label', size = 18)

# Confusion matrix
#cm = confusion_matrix(test_labels, rf_predictions)
cm = confusion_matrix(
    test_labels.argmax(axis=1), rf_predictions.argmax(axis=1))
plot_confusion_matrix(cm, classes = ["Bone","Lung","Liver"], title = 'Health Confusion Matrix')

plt.savefig('cm_BRCA.png')


# In[293]:


type(rf_predictions)


# In[294]:


type(test_labels)


# In[295]:


confusion_matrix(
    test_labels.argmax(axis=1), rf_predictions.argmax(axis=1))


# In[ ]:




