#!/usr/bin/env python
# coding: utf-8

# The objecitve of this notebook is to classify the TNM staging information for 8 of the cancer types. 
# We will add to the analysis the DE genes from each of the analyses completed in R

# We will extract features of the highest importance for each class.
# The catagorical classification will identify the TNM staging but will also illuminate features that
# define each class. 


# To comeplte the work, we will lean on the analysis completed on 4/7/2020
# We will construct and tune an RF to classify transcritomic data.

# We will move further and tune hypter paramters in a grid search for optimal mtrys and we will produce 
# ouput figures to visualize our work. 

# import libraries and packages necessary for analysis

# Load the Libraries

# import some packages thhat may be helpful for us

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

# define a random seed to replicate the results
RSEED = 50


# Read in the data set

features= pd.read_csv("/home/mskaro1/storage/PanCancerAnalysis/ML_2019/TNM_RF/TCGA-THCA_TNM.csv")

# Look at the shape of the data to construct the input layer. 
print('We have {} instances of data with {} variables'.format(*features.shape))



# dataframe.size 
size = features.size 

# dataframe.shape 
shape = features.shape 
  
# dataframe.ndim 
df_ndim = features.ndim 



size

shape

df_ndim


# In[329]:


# Labels are the values we want to predict, retrun and pop them off after
labels = np.array(features.pop('labels'))
# Saving feature names for later use
feature_list = list(features.columns)
# Convert to numpy array
features = np.array(features)


# In[330]:


type(labels)


# In[331]:


from numpy import *
labels=np.array(['NX' if x is np.nan else x for x in labels])

#i think this is turned into a list and not an np array which changes things, we may be able to deal with this
# earlier. 


# In[332]:


# split the data into train and test.

# 30% examples in test data
train, test, train_labels, test_labels = train_test_split(features,
                                         labels, 
                                         stratify = labels,
                                         test_size = 0.3, 
                                         random_state = RSEED)


# In[333]:


# Create the model with 1000 trees.

# Lets start with 1000 trees and see where it goes.


model = RandomForestClassifier(n_estimators=1000, 
                               random_state=RSEED, 
                               max_features = 'sqrt',
                               n_jobs=-1, verbose = 1)


# In[334]:


# Fit on training data. This will take some doing but seems viable. 

# from here down it is plug and play. I think fitting the data
# into memory may be an issue but we will see. 

model.fit(train, train_labels)


# In[335]:


# Training predictions (to demonstrate overfitting)
train_rf_predictions = model.predict(train)
train_rf_probs = model.predict_proba(train)[:, 1]

# Testing predictions (to determine performance)
rf_predictions = model.predict(test)
rf_probs = model.predict_proba(test)[:, 1]


# In[336]:


from sklearn.metrics import precision_score, recall_score, roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# Plot formatting
plt.style.use('fivethirtyeight')
plt.rcParams['font.size'] = 18


# In[337]:


'''def evaluate_model(predictions, probs, train_predictions, train_probs):
    """Compare machine learning model to baseline performance.
    Computes statistics and shows ROC curve."""
    
    # error in calculating baseline because of type error in the labels.
    # for loop won't work as at stands.
    # find solution to what needs to be done.
    baseline = {}
    
    baseline['recall'] = recall_score(test_labels, [1 for _ in range(len(test_labels))])
    baseline['precision'] = precision_score(test_labels, [1 for _ in range(len(test_labels))])
    baseline['roc'] = 0.5
    
    results = {}
    
    results['recall'] = recall_score(test_labels, predictions)
    results['precision'] = precision_score(test_labels, predictions)
    results['roc'] = roc_auc_score(test_labels, probs)
    
    train_results = {}
    train_results['recall'] = recall_score(train_labels, train_predictions)
    train_results['precision'] = precision_score(train_labels, train_predictions)
    train_results['roc'] = roc_auc_score(train_labels, train_probs)
    
    for metric in ['recall', 'precision', 'roc']:
        print(f'{metric.capitalize()} Baseline: {round(baseline[metric], 2)} Test: {round(results[metric], 2)} Train: {round(train_results[metric], 2)}')
    
    # Calculate false positive rates and true positive rates
    base_fpr, base_tpr, _ = roc_curve(test_labels, [1 for _ in range(len(test_labels))])
    model_fpr, model_tpr, _ = roc_curve(test_labels, probs)

    plt.figure(figsize = (8, 6))
    plt.rcParams['font.size'] = 16
    
    # Plot both curves
    plt.plot(base_fpr, base_tpr, 'b', label = 'baseline')
    plt.plot(model_fpr, model_tpr, 'r', label = 'model')
    plt.legend();
    plt.xlabel('False Positive Rate'); 
    plt.ylabel('True Positive Rate'); plt.title('ROC Curves');
    plt.show();
    plt.savefig('roc_auc_curve_BLCA.png');
    
    
    # save the results:
    
    
    
evaluate_model(rf_predictions, rf_probs, train_rf_predictions, train_rf_probs)
#plt.savefig('roc_auc_curve_BLCA.png');'''


# In[338]:


range(len(test_labels))


# In[339]:


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
cm = confusion_matrix(test_labels, rf_predictions)
plot_confusion_matrix(cm, classes = ['N0','NX','N1','N2','N3','M'])

plt.savefig('cm_THCA.png')


# In[340]:


print(type(rf_predictions))
print(type(rf_probs))
print(type(train_rf_predictions))
print(type(train_rf_probs))


# In[341]:


with open('TNM_output/THCA_rf_predictions.txt', 'w') as f:
    for item in rf_predictions:
        f.write("%s\n" % item)
with open('TNM_output/THCA_rf_probs.txt', 'w') as f:
    for item in rf_probs:
        f.write("%s\n" % item)
with open('TNM_output/THCA_train_rf_predictions.txt', 'w') as f:
    for item in train_rf_predictions:
        f.write("%s\n" % item)
with open('TNM_output/THCA_train_rf_probs.txt', 'w') as f:
    for item in train_rf_probs:
        f.write("%s\n" % item)


