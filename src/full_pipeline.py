#!/usr/bin/env python
# coding: utf-8

# In[175]:


# New approach lets attack building the shell.
import glob
import io
import math
import os
import re
import time
import warnings
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline
from imblearn.under_sampling import RandomUnderSampler, TomekLinks
from lightgbm import LGBMClassifier
from matplotlib import pyplot
from numpy import where
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_selection import RFE, SelectFromModel, SelectKBest, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import explained_variance_score
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import (
    LabelEncoder,
    MinMaxScaler,
    OneHotEncoder,
    OrdinalEncoder,
)

from numpy import loadtxt
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from sklearn.metrics import (
    classification_report,
    plot_confusion_matrix,
    plot_det_curve,
    plot_roc_curve,
)
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score


# In[176]:


class molecule_preprocessing:
    def __init__(self, path):
        self.path = path

        pass

    def CT(file):
        r = re.compile("TCGA-[a-zA-Z\d]{2,6}")
        m = r.search(file)
        if m:
            cancer_type = str(m.group(0))
        else:
            r = re.compile("TARGET-[a-zA-Z\d]{2,6}")
            m = r.search(file)
            cancer_type = str(m.group(0))
        cancer_type = cancer_type[5:]
        return cancer_type

    def read_selected(file, cancer_type):
        df = pd.read_table(file, delimiter=",")

        if "Unnamed: 0" in df.columns:
            df = df.drop("Unnamed: 0", axis=1)

        df = df[df.CT == cancer_type].fillna("NA")
        metastatic_sites = list(df.Metastatic_site)
        return metastatic_sites

    def read(file):

        df = pd.read_table(
            file,
            delimiter=",",
        ).fillna(0)

        if "barcode" in df.columns:
            df = df.drop("barcode", axis=1)
        if "__no_feature" in df.columns:
            df = df.drop("__no_feature", axis=1)
        if "Unnamed: 0" in df.columns:
            df = df.drop("Unnamed: 0", axis=1)
        if "__ambiguous" in df.columns:
            df = df.drop("__ambiguous", axis=1)
        if "__too_Low_aQual" in df.columns:
            df = df.drop("__too_low_aQual", axis=1)
        if "__not_aligned" in df.columns:
            df = df.drop("__not_aligned", axis=1)
        if "__alignment_not_unique" in df.columns:
            df = df.drop("__alignment_not_unique", axis=1)

        return df

    def cut(df, start, end, metastatic_site):

        X = df.iloc[:, start:end]

        # subsetted columns and the subsetted columns

        y = df[metastatic_site]

        # bind X block and y's
        X = pd.concat([X.reset_index(drop=True), y.reset_index(drop=True)], axis=1)

        return X

    def chunkIt(seq, num):
        """
        Split the feature sets into chunks
        """
        avg = len(seq) / float(num)
        out = []
        last = 0.0

        while last < len(seq):
            out.append(seq[int(last) : int(last + avg)])
            last += avg

        return out

    def split(df, metastatic_sites, site):
        for col in metastatic_sites:
            if col in df.columns:
                df[col] = df[col].replace(2, 1)
        y = df[site]

        for col in metastatic_sites:
            if col in df.columns:
                df = df.drop(col, axis=1)
        feature_list = df.columns
        return df, y, feature_list

    def synthetic_instances(X, y):
        counter = Counter(y_train)
        oversample = SMOTE()
        Xsm, ysm = oversample.fit_resample(X, y)
        Xsm = Xsm.astype(int)

        return Xsm, ysm


# In[177]:


class feature_selection:
    def __init__(self, X, y):
        self.X = X
        self.y = y
        pass

    def chi_selector(X, y, num_feats):
        feature_list = X.columns
        chi_selector = SelectKBest(chi2, k=num_feats)
        chi_selector.fit(X, y.ravel())
        chi_support = chi_selector.get_support()

        return chi_support, feature_list

    def rfR(X, y, num_feats):
        embeded_rf_selector = SelectFromModel(
            RandomForestRegressor(n_estimators=100, n_jobs=20), max_features=num_feats
        )
        embeded_rf_selector.fit(X, y.ravel())
        embeded_rf_support = embeded_rf_selector.get_support()
        # embeded_rf_feature = X.loc[:, embeded_rf_support].columns.tolist()

        return embeded_rf_support

    # chi_support = chi_selector(X, num_feats=5000)

    def logR(X, y, num_feats):
        rfe_selector = RFE(
            estimator=LogisticRegression(n_jobs=20),
            n_features_to_select=num_feats,
            step=40,
            verbose=0,
        )
        rfe_selector.fit(X, y.ravel())
        rfe_support = rfe_selector.get_support()

        return rfe_support

    # rfe_support = logR(X, num_feats=5000)

    def lassoR(X, y, num_feats):
        embeded_lr_selector = SelectFromModel(
            LogisticRegression(penalty="l2", n_jobs=20), max_features=num_feats
        )
        embeded_lr_selector.fit(X, y.ravel())
        embeded_lr_support = embeded_lr_selector.get_support()

        return embeded_lr_support

    # embeded_lr_support = lassoR(X, num_feats=5000)

    def rfC(X, y, num_feats):
        embeded_rf_selector = SelectFromModel(
            RandomForestClassifier(n_estimators=100, n_jobs=20), max_features=num_feats
        )
        embeded_rf_selector.fit(X, y.ravel())
        embeded_rf_support = embeded_rf_selector.get_support()

        return embeded_rf_support

    def cross_validate_feature_selection(
        feature_list,
        chi_support,
        rfe_support,
        embeded_lr_support,
        embeded_rfC_support,
        embeded_rfR_support,
    ):
        df = pd.DataFrame(
            {
                "Feature": feature_list,
                "Chi-2": chi_support,
                "RFE": rfe_support,
                "Logistics": embeded_lr_support,
                "RandomForestClassifier": embeded_rfC_support,
                "RandomForstRegression": embeded_rfR_support,
            }
        )
        # count the selected times for each feature
        df["Total"] = np.sum(df, axis=1)
        df = df.sort_values(["Total", "Feature"], ascending=False)
        df.index = range(1, len(df) + 1)

        return df

    def grade_features(X, y):
        chi_support, feature_list = feature_selection.chi_selector(X, y, num_feats=50)
        rfe_support = feature_selection.lassoR(X, y, num_feats=50)
        embeded_lr_support = feature_selection.logR(X, y, num_feats=50)
        embeded_rfC_support = feature_selection.rfC(X, y, num_feats=50)
        embeded_rfR_support = feature_selection.rfR(X, y, num_feats=50)

        CV = feature_selection.cross_validate_feature_selection(
            feature_list,
            chi_support,
            rfe_support,
            embeded_lr_support,
            embeded_rfC_support,
            embeded_rfR_support,
        )
        CV = CV[1:50]

        return CV


# In[76]:


mp = molecule_preprocessing(
    path="home/jovyan/storage/Clinical_annotation/Annotated_expression/*_2.csv"
)
file = "/home/jovyan/storage/Clinical_annotation/Annotated_expression/TCGA-ACC_annotated.csv"
cancer_type = molecule_preprocessing.CT(file=file)
selected_locs = molecule_preprocessing.read_selected(
    file="/home/jovyan/storage/Clinical_annotation/Sites_for_classiication.csv",
    cancer_type=cancer_type,
)
df = molecule_preprocessing.read(file=file)
dat = molecule_preprocessing.cut(df, start=0, end=500, metastatic_site=selected_locs[1])
X, y, feature_list = molecule_preprocessing.split(
    dat, selected_locs, site=selected_locs[1]
)
fs = feature_selection(X, y)
CV = feature_selection.grade_features(X, y)


# In[ ]:


import warnings

warnings.filterwarnings("ignore")


def grade_all_blocks(df, cancer_type, selected_locs, windows):
    for i in selected_locs:
        features = []
        for j in windows:
            # print(len(features))
            start = list(j)[0]
            end = list(j)[-1]
            if end <= df.shape[1]:
                sliced = molecule_preprocessing.cut(
                    df, start=start, end=end, metastatic_site=i
                )
                X, y, feature_list = molecule_preprocessing.split(
                    sliced, selected_locs, site=i
                )
                fs = feature_selection(X, y)
                CV = feature_selection.grade_features(X, y)
                # CV.to_csv("/home/jovyan/storage/Machine_Learning/Selected_features_binary_classifications/" + str("Organotropic_features")+ "_" + str(cancer_type) + "_" + str(i) + "_" + "ID_window_"  + str(start) + "_" + str(end) + ".csv")
                features.append(CV)

        dat = pd.concat(features, ignore_index=True)
        dat = dat.drop_duplicates(subset=["Feature"])
        dat = dat.sort_values(by="Total", ascending=False)
        dat.to_csv(
            "/home/jovyan/storage/Machine_Learning/Selected_features_binary_classifications/"
            + str("Organotropic_features")
            + "_"
            + str(cancer_type)
            + "_"
            + str(i)
            + ".csv"
        )
        pass


def grade_all_cancers(path):
    windows = molecule_preprocessing.chunkIt(seq=range(60483), num=100)
    for file in glob.glob(path):
        tic_all = time.perf_counter()
        print("Processing:" + str(file))
        cancer_type = molecule_preprocessing.CT(file=file)
        selected_locs = molecule_preprocessing.read_selected(
            file="/home/jovyan/storage/Clinical_annotation/Sites_for_classiication.csv",
            cancer_type=cancer_type,
        )
        df = molecule_preprocessing.read(file=file)
        grade_all_blocks(df, cancer_type, selected_locs, windows)
        print(str(file) + "Complete")
        toc_all = time.perf_counter()
        print(f"Read in features in {toc_all - tic_all:0.4f} seconds")
        pass


grade_all_cancers(
    path="/home/jovyan/storage/Clinical_annotation/Annotated_expression/*_2.csv"
)


# In[178]:


# Binary classification
class Organotropic_classification:
    def __init__(self, path):
        self.path = path
        
        pass
    
    
    def extract_info(file):
        # read the file name as a string
        
        # extract the cancer type
        r = re.compile("[A-Z\d]{3,8}")
        m = r.search(file)
        if m:
            CT = str(m.group(0))
        # extract the metastatic_site
        # re-name all of the files.
        
        r = re.compile(".*\_(.*)\.")
        m = r.search(file)
        
        if m:
            ms = str(m.group(1))
            if ms == "node":
                ms = "Lymph_node"
            if ms == "neck":
                ms = "Head_and_neck"
            if ms == "gland":
                ms = "Adrenal_gland"
            if ms == "tissue":
                ms = "Soft_tissue"
            if ms == "cavity":
                ms = "Oral_cavity"
            
        # extract the metastatic location
        df = pd.read_csv(file)
        features = df.Feature
        return features, CT, ms
    
    
    def subset(file, features, metastatic_site):
        # read in expression matrix
        df = pd.read_table(
            file,
            delimiter=",",
        ).fillna(0)

        if "barcode" in df.columns:
            df = df.drop("barcode", axis=1)
        if "__no_feature" in df.columns:
            df = df.drop("__no_feature", axis=1)
        if "Unnamed: 0" in df.columns:
            df = df.drop("Unnamed: 0", axis=1)
        if "__ambiguous" in df.columns:
            df = df.drop("__ambiguous", axis=1)
        if "__too_Low_aQual" in df.columns:
            df = df.drop("__too_low_aQual", axis=1)
        if "__not_aligned" in df.columns:
            df = df.drop("__not_aligned", axis=1)
        if "__alignment_not_unique" in df.columns:
            df = df.drop("__alignment_not_unique", axis=1)
        
        df[metastatic_site] = df[metastatic_site].replace(2, 1)
        y = df[metastatic_site]
        
        # select the the features
        df = df[df.columns & features]
        df = pd.concat([df.reset_index(drop=True), y.reset_index(drop=True)], axis=1)
        df = df.iloc[: , :-1]
        
        return df, y
    
    def d_split(X,y):
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.33, random_state=42)
        
        return X_train, X_test, y_train, y_test
    

    
    def synthetic_instances(X, y):
        '''
        make synthetic instances for classification this time
        '''
        print(sorted(Counter(y).items()))
        oversample = SMOTE()
        Xsm, ysm = oversample.fit_resample(X, y)
        print(sorted(Counter(ysm).items()))
        Xsm = Xsm.astype(int)

        return Xsm, ysm
    
    def classify(X_train, X_test, y_train, y_test):
        '''
        Input the data frames and put them into the classifier
        '''
        model = XGBClassifier()
        model.fit(X_train, y_train.ravel())
        y_pred = model.predict(X_test)
        predictions = [round(value) for value in y_pred]
        # evaluate predictions
        accuracy = accuracy_score(y_test, predictions)
        precision = precision_score(y_test, y_pred, average='weighted')
        recall = recall_score(y_test, y_pred, average='weighted')
        #f1 = f1_score(y_test, y_pred, average='weighted')
        print("Accuracy: %.2f%%" % (accuracy * 100.0))
        print("Precision: %.2f%%" % (precision * 100.0))
        print("Recall: %.2f%%" % (recall * 100.0))
        #print("F1 Score: %.2f%%" % (f1_score * 100.0))
        
        
        
        return model, X_test, y_test, y_pred, predictions
    
    
    def classify_2(X_train, X_test, y_train, y_test):
        '''
        Input the data frames and put them into the classifier
        '''
        
        model = RandomForestClassifier(n_estimators=1000, 
                               random_state=RSEED, 
                               max_features = 'sqrt',
                               n_jobs=-1, verbose = 1)
        
        model.fit(X_train, y_train.ravel())
        y_pred = model.predict(X_test)
        predictions = [round(value) for value in y_pred]
        # evaluate predictions
        accuracy = accuracy_score(y_test, predictions)
        precision = precision_score(y_test, y_pred, average='weighted')
        recall = recall_score(y_test, y_pred, average='weighted')
        #f1 = f1_score(y_test, y_pred, average='weighted')
        print("Accuracy: %.2f%%" % (accuracy * 100.0))
        print("Precision: %.2f%%" % (precision * 100.0))
        print("Recall: %.2f%%" % (recall * 100.0))
        #print("F1 Score: %.2f%%" % (f1_score * 100.0))
    
    
    
    def visualize(model, X_test, y_test, y_pred, predictions, ms, CT):
        from sklearn.metrics import (
            classification_report,
            plot_confusion_matrix,
            plot_det_curve,
            plot_roc_curve,
        )
        # may want to save the cfm here
        titles_options = [
        ("Confusion matrix, without normalization", None),
        ("Normalized confusion matrix", "true"),]
        for title, normalize in titles_options:
            disp = plot_confusion_matrix(
                model, X_test, y_test, cmap=plt.cm.Blues, normalize=normalize
            )
            disp.ax_.set_title(title)

            print(title)
            print(disp.confusion_matrix)
            plt.title(CT + " progresses to " + ms, fontweight="bold", fontsize=20)
            plt.xlabel("Metrics", fontweight="bold", fontsize=20)  # x-axis label with fontsize 15
            plt.ylabel("Classes", fontweight="bold", fontsize=20)  # y-axis label with fontsize 15
            plt.yticks(rotation=0)
            plt.savefig(CT + "_" + ms + "_" "reba_metrics_GBT_RNA_balanced.pdf")
        pass
    
    
    
    


# In[31]:


OC = Organotropic_classification(path = 
    "/home/jovyan/storage/Machine_Learning/Selected_features_binary_classifications/*.csv")


# In[ ]:


import warnings
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score

warnings.filterwarnings("ignore")

path = "/home/jovyan/storage/Machine_Learning/Selected_features_binary_classifications/*.csv"
ExpressionPath = "/home/jovyan/storage/Clinical_annotation/Completed_feature_selection/"

for file1 in glob.glob(path):
    features, CT, ms = Organotropic_classification.extract_info(file = file1)
    file2 = os.path.join(ExpressionPath, 'TCGA-'+ CT + '_annotated_2.csv')
    X, y = Organotropic_classification.subset(file = file2, 
                                          features=features,
                                          metastatic_site=ms)
    
    X_train, X_test, y_train, y_test = Organotropic_classification.d_split(X,y)
    #X_train,y_train = Organotropic_classification.synthetic_instances(X_train,y_train)
    #X_test,y_test = Organotropic_classification.synthetic_instances(X_test,y_test)
    
    if sorted(Counter(y).items())[1][1] >8:
        if X.shape[0] >50:
            if Counter(y_test)[1] >8:
                X_train,y_train = Organotropic_classification.synthetic_instances(X_train,y_train)
                X_test,y_test = Organotropic_classification.synthetic_instances(X_test,y_test)
                print(CT + "-" + ms)
                model, X_test, y_test, y_pred, predictions = Organotropic_classification.classify(X_train, X_test, y_train, y_test)
                Organotropic_classification.visualize(model, X_test, y_test, y_pred, predictions, ms, CT)
    
# DONE
