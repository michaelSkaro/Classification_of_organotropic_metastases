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
   
   # https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.chi2.html
   def chiSQ_selector(X, y, nFeatures):
       feature_list = X.columns
       chiSQ = SelectKBest(chi2, k= nFeatures)
       chiSQ.fit(X, y.ravel())
       chiSQ_support = chiSQ.get_support()

       return chiSQ-support, feature_list
   
   # https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html
   def rfR(X, y, nFeatures,nJobs):
       rfR_features = SelectFromModel(
           RandomForestRegressor(n_estimators=100, n_jobs = nJobs), max_features=nFeatures
       )
       rfR_features.fit(X, y.ravel())
       rfR_support = rfR_features.get_support()

       return rfR_support
   
   # https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFE.html
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
   
   
   # https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html
   def lassoR(X, y, nFeatures,nJobs):
       lR_selector = SelectFromModel(
           LogisticRegression(penalty="l2", n_jobs= nJobs), max_features=nFeatures
       )
       lR_selector.fit(X, y.ravel())
       lR_support = lR_selector.get_support()

       return lR_support
   
   
   # https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
   def rfC(X, y, nFeatures,nJobs):
       rfC_features = SelectFromModel(
           RandomForestClassifier(n_estimators=100, n_jobs=nJobs), max_features=nFeatures
       )
       rfC_features.fit(X, y.ravel())
       rfC_support = rfC_features.get_support()

       return rfC_support
   
   # cross validation 
   # https://towardsdatascience.com/the-5-feature-selection-algorithms-every-data-scientist-need-to-know-3a6b566efd2
   def cross_validate_feature_selection(
       feature_list,
       chiSQ_support,
       rfe_support,
       lR_support,
       rfC_support,
       rfR_support,
   ):
       df = pd.DataFrame(
           {
               "Feature": feature_list,
               "chi2": chiSQ_support,
               "RFE": rfe_support,
               "Logistics": lR_support,
               "RandomForestClassifier": rfC_support,
               "RandomForstRegression": rfR_support,
           }
       )
       
       
      df["Total"] = np.sum(df, axis=1)
       df = df.sort_values(["Total", "Feature"], ascending=False)
       df.index = range(1, len(df) + 1)

       return df

   def grade_features(X, y, nFeatures, n_jobs):
       chiSQ_support, feature_list = feature_selection.chiSQ_selector(X, y, nFeatures=nFeatures)
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

def main(cla_args=sys.argv[1:]):
   #--------------------------------#
   args = fs.args_parse(parse_cla(cla_args))
   X = args['x']
   y = args['y']
   nFeatures = args['f']
   nJobs = args['j']
   CV = feature_selection.grade_features(X = X, y = y, nFeatures = nFeatures, nJobs = nJobs)
   
   #--------------------------------#

if __name__ == '__main__':
    main()
