#!/usr/bin/env python
# coding: utf-8

# The objective of this notebook will be to read in the split data the runner has made
# learn on the split/annotated data and ake predictions using RnadomForest class in SciKitlearn

# Load the Libraries

# import some packages thhat may be helpful for us

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


# In[751]:


# define a random seed to replicate the results
RSEED = 50


# In[752]:


RSEED


# In[753]:


df= pd.read_csv("~/storage/PanCancerAnalysis/ML_2019/RF_input/TCGA-BRCA_DL_NDL.csv")


# In[707]:


print('We have {} pateints of data with {} variables'.format(*df.shape))

# dataframe.size 
size = df.size 
  
# dataframe.shape 
shape = df.shape 
  
# dataframe.ndim 
df_ndim = df.ndim 


# In[708]:


size


# In[709]:


shape


# In[710]:


df = df.rename(columns={'lables': 'labels'})
df = df.rename(columns={'labels': 'label'})


# In[711]:


labels = np.array(df.pop('label'))


# In[712]:


labels


# In[713]:


type(labels)


# In[714]:


from numpy import *


where_are_NaNs = isnan(labels)
labels[where_are_NaNs] = 0


# In[715]:


labels= labels.astype(int)


# In[716]:


type(labels)


# In[717]:


type(labels[1])


# In[719]:


labels.size


# In[722]:


# List of features for later use
feature_list = list(df.columns)


# In[724]:


# 30% examples in test data
train, test, train_labels, test_labels = train_test_split(df,
                                         labels, 
                                         stratify = labels,
                                         test_size = 0.3, 
                                         random_state = RSEED)


# In[725]:


# Imputation of missing values
train = train.fillna(train.mean())
test = test.fillna(test.mean())


# In[728]:


# Create the model with 1000 trees.

# Lets start with 1000 trees and see where it goes.


model = RandomForestClassifier(n_estimators=1000, 
                               random_state=RSEED, 
                               max_features = 'sqrt',
                               n_jobs=-1, verbose = 1)


# In[729]:


# Fit on training data. This will take some doing but seems viable. 

# from here down it is plug and play. I think fitting the data
# into memory may be an issue but we will see. 

model.fit(train, train_labels)


# In[731]:


# This is where we will customize the model. We will have to conduct a few tuning methods. 
# This will include changing the max feature, max depth, n estimators.
# I would like to also extract the importance features in each of the data. Kappa score might also be nice. 

# the first portion will be to find the best features for each class

# List of features for later use
# Get numerical feature importances
importances = list(model.feature_importances_)
# List of tuples with variable and importance
feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature_list, importances)]
# Sort the feature importances by most important first
feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
# Print out the feature and importances 
#[print('Variable: {:200} Importance: {}'.format(*pair)) for pair in feature_importances]


# In[638]:


n_nodes = []
max_depths = []

# Stats about the trees in random forest
for ind_tree in model.estimators_:
    n_nodes.append(ind_tree.tree_.node_count)
    max_depths.append(ind_tree.tree_.max_depth)
    
print(f'Average number of nodes {int(np.mean(n_nodes))}')
print(f'Average maximum depth {int(np.mean(max_depths))}')


# In[740]:


# list of x locations for plotting
x_values = list(range(len(importances)))

# List of features sorted from most to least important
sorted_importances = [importance[1] for importance in feature_importances]
sorted_features = [importance[0] for importance in feature_importances]
# Cumulative importances
cumulative_importances = np.cumsum(sorted_importances)

#print(np.where(cumulative_importances > 0.95))

# Find number of features for cumulative importance of 95%
# Add 1 because Python is zero-indexed
#print('Number of features for 95% importance:', np.where(cumulative_importances > 0.95)[0][0])


# In[741]:


print(sorted_importances[:10])


# In[742]:


# Training predictions (to demonstrate overfitting)
train_rf_predictions = model.predict(train)
train_rf_probs = model.predict_proba(train)[:, 1]

# Testing predictions (to determine performance)
rf_predictions = model.predict(test)
rf_probs = model.predict_proba(test)[:, 1]


# In[743]:


from sklearn.metrics import precision_score, recall_score, roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# Plot formatting
plt.style.use('fivethirtyeight')
plt.rcParams['font.size'] = 18


# In[744]:


def evaluate_model(predictions, probs, train_predictions, train_probs):
    """Compare machine learning model to baseline performance.
    Computes statistics and shows ROC curve."""
    
    baseline = {}
    
    baseline['recall'] = recall_score(test_labels, 
                                     [1 for _ in range(len(test_labels))])
    baseline['precision'] = precision_score(test_labels, 
                                      [1 for _ in range(len(test_labels))])
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
#plt.savefig('roc_auc_curve_BLCA.png');


# In[648]:


rf_predictions


# In[642]:


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
plot_confusion_matrix(cm, classes = ['Metastatic Cancer', 'Non-metastatic Cancer'],
                      title = 'Health Confusion Matrix')

plt.savefig('cm_BRCA.png')




