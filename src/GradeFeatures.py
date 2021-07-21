#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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


# In[20]:


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
    
    def cut(df, start,end,metastatic_sites):
        
        X = df.iloc[:, start:end]

        # subsetted columns and the subsetted columns

        y = df[metastatic_sites]

        sites = list(y.columns)

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
        y = df[site]
        for col in metastatic_sites:
            if col in df.columns:
                df[col] = df[col].replace(2, 1)
        for col in metastatic_sites:
            if col in df.columns:
                df = df.drop(col, axis=1)
        X = df
        feature_list = X.columns
        return X, y, feature_list

    def synthetic_instances(X, y):
        counter = Counter(y_train)
        oversample = SMOTE()
        Xsm, ysm = oversample.fit_resample(X, y)
        Xsm = Xsm.astype(int)

        return Xsm, ysm


# In[21]:


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


# In[25]:


mp = molecule_preprocessing(
    path="home/jovyan/storage/Clinical_annotation/Annotated_expression/*_2.csv"
)
file = "/home/jovyan/storage/Clinical_annotation/Annotated_expression/TCGA-LUAD_annotated.csv"
cancer_type = molecule_preprocessing.CT(file=file)
selected_locs = molecule_preprocessing.read_selected(
    file="/home/jovyan/storage/Clinical_annotation/Sites_for_classiication.csv",
    cancer_type=cancer_type,
)
df = molecule_preprocessing.read(
    file=file)
dat = molecule_preprocessing.cut(df, start=0,end=500,metastatic_sites=selected_locs)
X, y, feature_list = molecule_preprocessing.split(
    dat, selected_locs, site=selected_locs[1]
)
fs = feature_selection(X, y)
CV = feature_selection.grade_features(X, y)


# In[ ]:


import warnings
warnings.filterwarnings('ignore')

def grade_all_blocks(df,cancer_type,selected_locs,windows):
    for i in selected_locs:
        features = []
        for j in windows:
            start = list(j)[0]
            end = list(j)[-1]
            # this is what is slowing everything down.
            sliced = molecule_preprocessing.cut(df, start=start,end=end,metastatic_sites=selected_locs)
            X, y, feature_list = molecule_preprocessing.split(
                sliced, selected_locs, site=i
            )
            fs = feature_selection(X, y)
            CV = feature_selection.grade_features(X, y)
            # CV.to_csv("/home/jovyan/storage/Machine_Learning/Selected_features_binary_classifications/" + str("Organotropic_features")+ "_" + str(cancer_type) + "_" + str(i) + "_" + "ID_window_"  + str(start) + "_" + str(end) + ".csv")
            features.append(CV)

        dat = pd.concat(features, ignore_index=True)
        dat = dat.drop_duplicates(subset=["Feature"])  # ****
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
    windows = molecule_preprocessing.chunkIt(seq=range(60485), num=100)
    for file in glob.glob(path):
        cancer_type = molecule_preprocessing.CT(file=file)
        selected_locs = molecule_preprocessing.read_selected(
            file="/home/jovyan/storage/Clinical_annotation/Sites_for_classiication.csv",
            cancer_type=cancer_type,
        )
        df = molecule_preprocessing.read(
        file=file)
        tic_all = time.perf_counter()
        print("Processing:" + str(file))
        grade_all_blocks(df,cancer_type,selected_locs,windows)
        print(str(file) + "Complete")
        toc_all = time.perf_counter()
        print(f"Read in features in {toc_all - tic_all:0.4f} seconds")
        pass


grade_all_cancers(
    path="/home/jovyan/storage/Clinical_annotation/Annotated_expression/*_2.csv"
)


# In[ ]:




