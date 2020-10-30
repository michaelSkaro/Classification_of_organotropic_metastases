###############################################################################
# 10/28/20

'''
Objective: Figure quality ROC and PRC curves for MOT project 
Input:  Raw data, classifications, model outputs for plotting
Output: Visualization + code to integrate into the pipeline
'''
###############################################################################
# my major imports 
from sklearn.model_selection import train_test_split
import numpy as np
rom sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import plot_precision_recall_curve
import matplotlib.pyplot as plt
###############################################################################

###############################################################################
#//Input model here:
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

# import whatever DF you need too, must have labels as last column, features across the topp and 
# the instances down the side, get rid of row names if you need too

X= pd.read_csv("~/storage/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/TCGA-BRCA_metastatic_data_RNAseq.csv")
X.drop(X.columns[0], axis=1)
# Look at the shape of the data to construct the input layer. 
print('We have {} instances of data with {} variables'.format(*X.shape))


y = X[["Bone","Lung","Liver"]].to_numpy()
n_classes = y.shape[1]
patients = np.array(X.pop('barcode'))
# remove the features labels
X = X.drop(["Bone","Lung","Liver"], axis = 1) 

# Labels are the values we want to predict, retrun and pop them off after
# List of features for later use
feature_list = list(X.columns)
# Convert to numpy array
features = np.array(X)


X_train, X_test, y_train, y_test = train_test_split(X,
                                         y, 
                                         stratify = labels,
                                         test_size = 0.3, 
                                         random_state = RSEED)
# choose the model:RF classification
import random
RSEED = random.randint()
model = RandomForestClassifier(n_estimators=100, 
                               random_state=RSEED, 
                               max_features = 'sqrt', 
                               n_jobs=-1, verbose = 1)
# Fit training data
# this may change once marcus uses the IG algo, gpotta discuss today/tomorrow?, Also are we wrapping in Onevs rest or what?
model.fit(X_train, y_train) # or do we wrap before fit, or is it a different fit call? look at sklearn.documentation

###############################################################################

# Compute the average precision score

from sklearn.metrics import average_precision_score
average_precision = average_precision_score(y_test, y_score)

print('Average precision-recall score: {0:0.2f}'.format(
      average_precision))
      
      
###############################################################################
# Plot the Precision-Recall curve: Binary 


disp = plot_precision_recall_curve(classifier, X_test, y_test)
disp.ax_.set_title('Binary Precision-Recall curve: '
                   'AP={0:0.2f}'.format(average_precision)) # input for each stage of work
                   
###############################################################################
# The average precision score in multi-label settings

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score

# For each class
precision = dict()
recall = dict()
average_precision = dict()
for i in range(n_classes):
    precision[i], recall[i], _ = precision_recall_curve(y_test[:, i],
                                                        y_score[:, i])
    average_precision[i] = average_precision_score(y_test[:, i], y_score[:, i])

# A "micro-average": quantifying score on all classes jointly
precision["micro"], recall["micro"], _ = precision_recall_curve(y_test.ravel(),
    y_score.ravel())
average_precision["micro"] = average_precision_score(y_test, y_score,
                                                     average="micro")
print('Average precision score, micro-averaged over all classes: {0:0.2f}'
      .format(average_precision["micro"]))
      
###############################################################################
# Plot the micro-averaged Precision-Recall curve

plt.figure()
plt.step(recall['micro'], precision['micro'], where='post')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title(
    'Average precision score, micro-averaged over all classes: AP={0:0.2f}'
    .format(average_precision["micro"]))      


###############################################################################
# Plot Precision-Recall curve for each class and iso-f1 curves

#
from itertools import cycle
# setup plot details
colors = cycle(["red","black","blue","green","darkorange","teal","turquoise","navy"]) 
# we need to check how many classes we need to plot here though, I think BLCA has 7 
# Can change for the different comparisons but i think this is probably fine. 

plt.figure(figsize=(7, 8))
f_scores = np.linspace(0.2, 0.8, num=4)
lines = []
labels = []
for f_score in f_scores:
    x = np.linspace(0.01, 1)
    y = f_score * x / (2 * x - f_score)
    l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
    plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

lines.append(l)
labels.append('iso-f1 curves')
l, = plt.plot(recall["micro"], precision["micro"], color='gold', lw=2)
lines.append(l)
labels.append('micro-average Precision-recall (area = {0:0.2f})'
              ''.format(average_precision["micro"]))

for i, color in zip(range(n_classes), colors):
    l, = plt.plot(recall[i], precision[i], color=color, lw=2)
    lines.append(l)
    labels.append('Precision-recall for class {0} (area = {1:0.2f})'
                  ''.format(i, average_precision[i]))

fig = plt.gcf()
fig.subplots_adjust(bottom=0.25)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Extension of Precision-Recall curve to multi-class')
plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=14))


plt.show()     

# 10/29/20

###############################################################################
#Run model code above, to continue***
###############################################################################




# Major inports i need now, probs need more later... 
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.metrics import average_precision_score

from sklearn.ensemble import RandomForestClassifier
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier

from numpy import interp # numpy.interp()**** scipy one is depricated...

from sklearn.metrics import roc_auc_score # for the area under curve scoring

# Compute ROC curve and ROC area for each class
from sklearn.preprocessing import label_binarize
lab = label_binarize(labels, classes=unique(labels))
n_classes = shape(lab[1])

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

###############################################################################
#Plot of a ROC curve for  just one class
###############################################################################

plt.figure()
lw = 2
plt.plot(fpr[2], tpr[2], color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()


##############################################################################
# Plot ROC curves for the multilabel output

# Compute macro-average ROC curve and ROC area
# aggregate all false positive rates
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

# interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += interp(all_fpr, fpr[i], tpr[i])

# average it and compute AUC
mean_tpr /= n_classes

fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# Plot all ROC curves
plt.figure()
plt.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]),
         color='deeppink', linestyle=':', linewidth=4)

plt.plot(fpr["macro"], tpr["macro"],
         label='macro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle=':', linewidth=4)

colors = cycle(["red","black","blue","green","darkorange","teal","turquoise","navy"]) 
#same number colors for the PRCs. Although now that I am getting into this, I think the binaries will be fine?

for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))

plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic to multi-class')
plt.legend(loc="lower right")
plt.show()

# Done
