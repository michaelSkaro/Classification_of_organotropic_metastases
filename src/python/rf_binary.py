'''
Random Forest implementation with
code for hyper-parameter tuning.
'''


import os
import re
import sys
import csv
import argparse


import numpy as np
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, \
                            confusion_matrix

def parse_cli(args):
    '''
    Parses command line arguments.

    Parameters
    ----------
    args: list
        Command line options followed by arguments.

    Returns
    -------
    parsed_args: dict
        Key-value pairs of options and arguments.
    '''

    parser = argparse.ArgumentParser(
        description='Performs binary classification using a Random Forest' \
                    + ' Classifier.'
    )
    parser.add_argument('--train', required=True, 
                        help='Path to the training dataset (.csv).')
    parser.add_argument('--test', required=True, 
                        help='Path to the testing dataset (.csv).')
    parser.add_argument('-o', '--outdir', required=False, default=None, 
                        help='Path to the directory to store results.')

    return vars(parser.parse_args(args))
    

def classify(x_train, y_train, x_test, y_test, attrs, output_file_path=None):

    clf = RandomForestClassifier(oob_score=True)
    clf.fit(x_train, y_train)
    preds = clf.predict(x_test)

    report = classification_report(y_test, preds, 
                                   labels=['Positive', 'Negative'])
    cf_matrix = confusion_matrix(y_test, preds, 
                                         labels=['Positive', 'Negative'])
    
    if output_file_path is None:
        print('Model Training OOB Error:', 1 - clf.oob_score_)
        print('\nClassification Report:\n', report) 
        print('Confusion Matrix:\n', cf_matrix)

        #Feature Importance
        important_features = np.column_stack([
            attrs[:-1], 
            clf.feature_importances_
        ])
        inds = np.argsort(important_features[:, -1])
        inds = np.flip(inds)
        important_features = important_features[inds]
        print('\nImportant Features:\n', important_features[:5])

        '''
        for i, estimator in enumerate(clf.estimators_):
            print('Decision Tree', i,  ' Node Depth:', estimator.tree_.max_depth,
                  ' Num of Leaves:', estimator.tree_.n_leaves)
        '''
    else:
        with open(output_file_path, 'w+') as fp:
            fp.write('Classification Report:\n')
            fp.write(report)
            fp.write('\n Confusion Matrix:\n')
            fp.write(np.array2string(cf_matrix))

def plot_training_oob(x_train, y_train):
    '''
    Plots the Out-of-Bag score of the Random Forest.

    Results: n_estimators: 500, max_depth: 4, max_leaves: 8
    '''
    
    errors = []
    #num_of_estimators = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
    #depths = np.arange(1, 11) 
    n_leaves = np.arange(2, 15) 
    for max_leaf_nodes in tqdm(n_leaves):
        clf = RandomForestClassifier(n_estimators=100, oob_score=True,
                                     max_depth=None, max_leaf_nodes=max_leaf_nodes)
        clf.fit(x_train, y_train)
        errors.append(1 - clf.oob_score_)

    sns.set()
    fig = plt.figure()
    plt.suptitle('Random Forest Training Loss')
    ax = fig.add_subplot()
    ax.set_ylabel('OOB Error')
    ax.set_xlabel('Max Leaf Nodes')
    sns.lineplot(x=n_leaves, y=errors, ax=ax)
    plt.show()

def main():

    args = parse_cli(sys.argv[1:])
    train_path = args['train']
    test_path = args['test']
    outdir_path = args['outdir']

    with open(train_path) as fp:
        reader = csv.reader(fp)
        train_data = list(reader)

    with open(test_path) as fp:
        reader = csv.reader(fp)
        test_data = list(reader)

    train_data = np.array(train_data)
    test_data = np.array(test_data)
    x_train, y_train = train_data[1:, :-1].astype(np.float), train_data[1:, -1]
    x_test, y_test = test_data[1:, :-1].astype(np.float), test_data[1:, -1]

    if outdir_path is None:
        classify(x_train, y_train, x_test, y_test, test_data[0])
    else:
        output_file_name = os.path.split(train_path)[1]
        output_file_name = re.sub('train', 'rf_result', output_file_name)[:-4] + '.txt'
        output_file_path = os.path.join(outdir_path, output_file_name)
        classify(x_train, y_train, x_test, y_test, test_data[0], output_file_path)

    #plot_training_oob(x_train, y_train)


if __name__ == '__main__':
    main()
