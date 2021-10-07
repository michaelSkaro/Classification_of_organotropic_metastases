#!/usr/local/bin/python3

# Make Shap.py

# imports
import sys
import pandas as pd
import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
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
import shap
#
#


# Load data, we will make this an arguement at the end

def extract_info(file, target):
    df = pd.read_csv(file, delimiter = "\t") # Load the data, make the command line arguements for the model
    
    # get the x and y variables back from the data to put into basic model splitting for shit model
    y = df[target]
    
    X = df.drop(target, axis=1)
    
    return X,y

def d_split(X,y):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.33, random_state=42)

    return X_train, X_test, y_train, y_test

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

    return model, X_test, y_test, y_pred, predictions, accuracy, precision, recall

# shap this bitch

def shapperBAR(model, X_test, y_test, y_pred, predictions):
    shap_values = shap.TreeExplainer(model).shap_values(X_train)
    f = plt.figure()
    shap.summary_plot(shap_values, X_train, plot_type="bar")
    f.savefig(col + "_bar.png", bbox_inches='tight', dpi=600)
    return shap_values

# not complete

def mayTheShapperForcebeWithYou(model, X_test, y_test, y_pred, predictions, j):
    shap.initjs()
    # Get the predictions and put them with the test data.
    X_output = X_test.copy()
    X_output.loc[:,'predict'] = np.round(model.predict(X_output),2)

    # Randomly pick some observations
    random_picks = np.arange(1,330,50) # Every 50 rows
    S = X_output.iloc[random_picks]
    explainerModel = shap.TreeExplainer(model)
    shap_values_Model = explainerModel.shap_values(S)
    p = shap.force_plot(explainerModel.expected_value, shap_values_Model[j], S.iloc[[j]])
    return(p)

def makeDependence(X_train, shap_values):
	
	for col in X_train.columns:
		
		for i in len(shap.values):
			f = plt.figure()
			shap.dependence_plot(col, shap_values[1], X_train)
			f.savefig(col + "_dependence.png", bbox_inches='tight', dpi=600)

	pass




def main():
	args = sys.argv[1:]
    df = extract_info(file = args[0], target = args[1])
	X_train, X_test, y_train, y_test = d_split(X,y)
	model, X_test, y_test, y_pred, predictions, accuracy, precision, recall = classify(X_train, X_test, y_train, y_test)
	shap_values = shapperBAR(model, X_test, y_test, y_pred, predictions)
	makeDependence(X_train, shap_values)

if __name__ == "__main__":
    main()

