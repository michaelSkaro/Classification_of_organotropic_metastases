#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


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


# In[3]:


import os
import re
import sys
import argparse
import subprocess
import math

import fs.args_parse

class feature_selection:
   def __init__(self, X, y, nFeatures,nJobs):
       self.X = X
       self.y = y
       self.nFeatures = nFeatures
       self.nJobs = nJobs
       pass
   # adapted from here: https://www.kaggle.com/mlwhiz/feature-selection-using-football-data
   # https://towardsdatascience.com/the-5-feature-selection-algorithms-every-data-scientist-need-to-know-3a6b566efd2
   
   def chiSQ-selector(X, y, nFeatures):
       feature_list = X.columns
       chiSQ-selector = SelectKBest(chi2, k= nFeatures)
       chiSQ-selector.fit(X, y.ravel())
       chiSQ-support = chiSQ-selector.get_support()

       return chiSQ-support, feature_list
   # adapted from here: https://www.kaggle.com/mlwhiz/feature-selection-using-football-data
   # https://towardsdatascience.com/the-5-feature-selection-algorithms-every-data-scientist-need-to-know-3a6b566efd2
   
   def rfR(X, y, nFeatures,nJobs):
       rfR_features = SelectFromModel(
           RandomForestRegressor(n_estimators=100, n_jobs = nJobs), max_features=nFeatures
       )
       rfR_features.fit(X, y.ravel())
       rfR_support = rfR_features.get_support()

       return rfR_support
   # adapted from here: https://www.kaggle.com/mlwhiz/feature-selection-using-football-data
   # https://towardsdatascience.com/the-5-feature-selection-algorithms-every-data-scientist-need-to-know-3a6b566efd2
   
   def recursiveFeatureSelection(X, y, nFeatures,nJobs):
       rfe_selector = RFE(
           estimator=LogisticRegression(n_jobs = nJobs),
           n_features_to_select=nFeatures,
           step= math.celi(len(X.columns)/10),
           verbose=0,
       )
       rfe_selector.fit(X, y.ravel())
       rfe_support = rfe_selector.get_support()

       return rfe_support
   
   # adapted from here: https://www.kaggle.com/mlwhiz/feature-selection-using-football-data
   # https://towardsdatascience.com/the-5-feature-selection-algorithms-every-data-scientist-need-to-know-3a6b566efd2
   
   def lassoR(X, y, nFeatures,nJobs):
       lR_selector = SelectFromModel(
           LogisticRegression(penalty="l2", n_jobs= nJobs), max_features=nFeatures
       )
       lR_selector.fit(X, y.ravel())
       lR_support = lR_selector.get_support()

       return lR_support
   
   # adapted from here: https://www.kaggle.com/mlwhiz/feature-selection-using-football-data
   # https://towardsdatascience.com/the-5-feature-selection-algorithms-every-data-scientist-need-to-know-3a6b566efd2
   
   def rfC(X, y, nFeatures,nJobs):
       rfC_features = SelectFromModel(
           RandomForestClassifier(n_estimators=100, n_jobs=nJobs), max_features=nFeatures
       )
       rfC_features.fit(X, y.ravel())
       rfC_support = rfC_features.get_support()

       return rfC_support
   
   def cross_validate_feature_selection(
       feature_list,
       chiSQ-support,
       rfe_support,
       lR_support,
       rfC_support,
       rfR_support,
   ):
       df = pd.DataFrame(
           {
               "Feature": feature_list,
               "chi2": chiSQ-support,
               "RFE": rfe_support,
               "Logistics": lR_support,
               "RandomForestClassifier": rfC_support,
               "RandomForstRegression": rfR_support,
           }
       )
       # adapted from here: https://www.kaggle.com/mlwhiz/feature-selection-using-football-data
       # https://towardsdatascience.com/the-5-feature-selection-algorithms-every-data-scientist-need-to-know-3a6b566efd2
       
      df["Total"] = np.sum(df, axis=1)
       df = df.sort_values(["Total", "Feature"], ascending=False)
       df.index = range(1, len(df) + 1)

       return df

   def grade_features(X, y, nFeatures, n_jobs):
       chiSQ-support, feature_list = feature_selection.chiSQ-selector(X, y, nFeatures=nFeatures)
       rfe_support = feature_selection.recursiveFeatureSelection(X, y, nFeatures=nFeatures, nJobs = nJobs)
       lR_support = feature_selection.lassoR(X, y, nFeatures=nFeatures, nJobs = nJobs)
       rfC_support = feature_selection.rfC(X, y, nFeatures=nFeatures, nJobs = nJobs)
       rfR_support = feature_selection.rfR(X, y, nFeatures=nFeatures, nJobs = nJobs)

       CV = feature_selection.cross_validate_feature_selection(
           feature_list,
           chiSQ-support,
           rfe_support,
           lR_support,
           rfC_support,
           rfR_support,
       )
       CV = CV[1:nFeatures]

       return CV

# In[4]:


mp = molecule_preprocessing(
    path="home/jovyan/storage/Clinical_annotation/Annotated_expression/*_2.csv"
)
file = "/home/jovyan/storage/Clinical_annotation/Annotated_expression/TCGA-ACC_annotated.csv"
cancer_type = molecule_preprocessing.CT(file=file)
selected_locs = molecule_preprocessing.read_selected(
    file="/home/jovyan/storage/Clinical_annotation/Sites_for_classiication.csv",
    cancer_type=cancer_type,
)
df = molecule_preprocessing.read(
    file=file)
dat = molecule_preprocessing.cut(df, start=0,end=500,metastatic_site=selected_locs[1])
X, y, feature_list = molecule_preprocessing.split(
    dat, selected_locs, site=selected_locs[1]
)
fs = feature_selection(X, y, 50, 20)
CV = feature_selection.grade_features(X, y, 50, 20)


# In[ ]:


import warnings

warnings.filterwarnings("ignore")


def grade_all_blocks(df, cancer_type, selected_locs, windows):
    for i in selected_locs:
        features = []
        for j in windows:
            #print(len(features))
            start = list(j)[0]
            end = list(j)[-1]
            if end<=df.shape[1]:
                sliced = molecule_preprocessing.cut(
                    df, start=start, end=end, metastatic_site=i
                )
                X, y, feature_list = molecule_preprocessing.split(
                    sliced, selected_locs, site=i
                )
                fs = feature_selection(X, y)
                CV = feature_selection.grade_features(X, y, 50, 20)
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

