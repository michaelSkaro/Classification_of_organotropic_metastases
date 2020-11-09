'''
Multi-label classifier of tumor metastasis locactions.

Dataset: Metastatic Loci One-hot encoded.
'''

import os
import csv
import sys
import argparse

import numpy as np
from sklearn.model_selection import RepeatedKFold
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, \
                            recall_score, f1_score

def load_data(file_path):
    '''
    Loads in data from a csv file and returns a NumPy 
    array.

    Parameters
    ----------
    file_path: string
        Path to the csv file.

    Returns
    -------
    data: NumPy array (n_rows, n_cols)
    '''

    data = []
    with open(file_path) as fp:
        reader = csv.reader(fp)
        data = list(reader)

    return np.array(data)

def feature_label_split(data):
    '''
    Splits the dataset into two sub-matrices: 
    a feature matrix and a label matrix.

    Paramters
    ---------
    data: NumPy Array (n_rows, n_cols):
        Dataset to be split.

    Returns
    -------
    x: NumPy Array (n_rows, n_cols - n_labels)
        Feature matrix.
    y: NumPy Array (n_rows, n_labels)
        Label matrix.
    feature_names: NumPy Array (n_features,)
        The name of each feature.
    labels: NumPy Array (n_labels,)
        Names of the labels.
    '''

    label_indices = [-1 if 'ENSG' in x  or 'embedding' in x else i 
                     for i, x in enumerate(data[0])]
    labels_starting_index = np.min([x for x in label_indices if x > 1])
    num_of_labels = len(data[0]) - labels_starting_index
    x = data[1:, 2:labels_starting_index].astype(np.float)
    y = data[1:, labels_starting_index:].astype(np.int)
    feature_names = data[0, 2:labels_starting_index]
    labels = data[0, labels_starting_index:]

    return x, y, feature_names, labels

def feature_importance(feature_names, model):
    '''
    Returns the most important features used by
    the classifier, in sorted order.

    Parameters
    ----------
    feature_names: NumPy array (n_features,)
        The names of all features in the dataset.
    model: A scikit-learn Estimator
        An estimator that contains a feature_importance_
        attribute.

    Returns
    -------
    sorted_featrues: NumPy Array (k, 2) 
        An array of the sorted features and their 
        importance value. 0 < k <= n_features.
    '''

    important_features = []
    for i, x in enumerate(model.feature_importances_): 
        if x > 0:
            important_features.append([feature_names[i], x])

    important_features = np.array(important_features)
    idx = np.argsort(important_features[:,-1])
    sorted_features = important_features[idx]
    
    return sorted_features

def write_feature_importance(model, labels, feature_names, outfile_path):
    '''
    Writes the important features of each estimator into a
    csv file.


    Parameters
    ----------
    model: A scikit-learn Estimator
        An estimator that contains a feature_importance_
        attribute.
    labels: NumPy Array (n_labels,) 
        Labels of the classses in the dataset.
    feature_names: NumPy Array (n_features,)
        The name of each feature.
    outfile_path: string
        Path to save the csv file.

    Returns
    -------
    None
    '''

    important_features = []
    for i, estimator in enumerate(model.estimators_):
        important_features.append(list())
        important_features[-1].append(labels[i])
        for row in feature_importance(feature_names, estimator):
            important_features[-1].append(row[0])

        important_features.append(list())
        important_features[-1].append('Importance Score')
        for row in feature_importance(feature_names, estimator):
            important_features[-1].append(row[1])

    longest_length = max([len(x) for x in important_features])
    with open(outfile_path, 'w') as fp:
        writer = csv.writer(fp, delimiter=',')
        for j in range(longest_length):
            output = []
            for i in range(len(important_features)):    
                if not j >= len(important_features[i]):
                    output.append(important_features[i][j])   

            writer.writerow(output)

def parse_cli(args):
    '''
    Parsers the command line arguments.
    
    Parameters
    ----------
    args: list
        List of argument-value pairs.

    Returns
    -------
    parsed_args: dict
        Parsed into a dictionary argument-value pairs.
    '''
    
    parser = argparse.ArgumentParser(
        description='Multi-label classifier of tumor metastasis locactions.'
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help='Directory path of csv dataset files.')
    parser.add_argument('-o', '--output_dir', required=False,
                        help='Path to a directory where the important ' + 
                             'features used in the models will be stored ' + 
                             'in csv files.')

    return vars(parser.parse_args(args))
    
def main():
    args = parse_cli(sys.argv[1:])
    file_dir_path = args['input_dir']
    for file_name in os.listdir(file_dir_path):
        if '.csv' in file_name and 'KIRP' not in file_name:
            print(file_name)
            file_path = os.path.join(file_dir_path, file_name)
            data = load_data(file_path)
            x, y, feature_names, labels = feature_label_split(data)
            num_of_instances = len(x)
            num_of_labels = len(labels)

            print('Number of instances:', num_of_instances)
            print('Number of labels:', num_of_labels)
            
            multi_label_accuracy = []
            class_wise_accuracy = [list() for i in range(num_of_labels)]
            class_wise_precision = [list() for i in range(num_of_labels)]
            class_wise_recall = [list() for i in range(num_of_labels)]
            class_wise_f1 = [list() for i in range(num_of_labels)]
            rkf = RepeatedKFold(n_splits=10, n_repeats=10, random_state=None)
     
            for train_indices, test_indices in rkf.split(x):
                x_train, y_train = x[train_indices], y[train_indices]
                x_test, y_test = x[test_indices], y[test_indices]
                model = OneVsRestClassifier(RandomForestClassifier(), n_jobs=-1)
                model.fit(x_train, y_train)
                preds = model.predict(x_test)
                multi_label_accuracy.append(accuracy_score(y_test, preds))

                for i in range(num_of_labels):
                    class_wise_accuracy[i].append(
                        accuracy_score(y_test[:, i], preds[:, i])
                    )

                    class_wise_precision[i].append(
                        precision_score(y_test[:, i], preds[:, i], 
                                        zero_division=0)
                    )

                    class_wise_recall[i].append(
                        recall_score(y_test[:, i], preds[:, i],
                                     zero_division=0)
                    )

                    class_wise_f1[i].append(
                        f1_score(y_test[:, i], preds[:, i],
                                 zero_division=0)
                    )
            
            mean_multi_label_accuracy = np.mean(multi_label_accuracy)
            mean_class_wise_accuracy = [np.mean(class_wise_accuracy[i]) 
                                        for i in range(num_of_labels)]
            mean_class_wise_precision = [np.mean(class_wise_precision[i]) 
                                         for i in range(num_of_labels)]
            mean_class_wise_recall = [np.mean(class_wise_recall[i]) 
                                      for i in range(num_of_labels)]
            mean_class_wise_f1 = [np.mean(class_wise_f1[i]) 
                                  for i in range(num_of_labels)]

            print('Average Multi-label Accuracy:', mean_multi_label_accuracy)
            print('Average Class-wise Metrics:')
            for i in range(num_of_labels):
                print('\t', labels[i] + ' - ', 
                      'Accuracy: ', mean_class_wise_accuracy[i],
                      'Precision: ', mean_class_wise_precision[i],
                      'Recall: ', mean_class_wise_recall[i],
                      'F1: ', mean_class_wise_f1[i],
                )

            print() #Separates the metrics of different datasets 
            
            if args['output_dir'] is not None:
                outfile_path = file_name[:-4] + '_important_features.csv'
                outfile_path = os.path.join(args['output_dir'], outfile_path)
                write_feature_importance(model, labels, feature_names,
                                         outfile_path)

if __name__ == '__main__':
    main()
