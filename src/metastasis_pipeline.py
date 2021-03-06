'''
A supervisory script to run all data through the
metastatic loci classificaiton pipeline.

Pipeline Steps:
    1. Create binary files.
    2. Oversample using SMOTE.
    3. Perform feature selection using Gain Ratio rankings.
    4. Classify (DT, RF, NN, SVM)
'''


import os
import re
import sys
import argparse
import subprocess

from tqdm import tqdm

import mot.smote
import mot.rf_binary
import mot.create_binary_datasets
import mot.gain_ratio_feature_selection


def parse_cli(args):
    '''
    Parses command line arguments.

    Parameters
    ----------
    args: list
        List of command line options followed by arguments.
    
    Returns
    -------
    args: dict
        Parsed arguments in key-value pairs.
    '''
    
    parser = argparse.ArgumentParser(
        description='A supervisory script to run all data through the'
                    + ' metastatic loci classificaiton pipeline.'
    )
    required = parser.add_argument_group('required args')
    required.add_argument('-i', '--input', required=True,
                          help='Path to input directory containing datasets.')
    required.add_argument('-o', '--output', required=True,
                          help='Path to output directory to save results.')
    required.add_argument('-w', '--weka_jar', required=True,
                          help='Path to Weka jar file.')
    required.add_argument('-c', '--classes_dir', required=True,
                          help='Path to save java bytecode.')
    required.add_argument('-j', '--java_src', required=True,
                          help='Path to the GainRatio java file.')

    return vars(parser.parse_args(args))

def create_directories(output_dir_path):
    '''
    Creates directories to store pipeline results.

    Parameters
    ----------
    output_dir_path: string
        Path to output directory to make new directories.

    Returns
    -------
    paths: tuple of strings
        Paths of all the directories created. 
    '''

    binary_dir_path = os.path.join(output_dir_path, 'binary-datasets')
    oversampled_dir_path = os.path.join(output_dir_path, 
                                        'oversampled-datasets')
    important_features_dir_path = os.path.join(output_dir_path,
                                               'important-features')
    feature_selected_dir_path = os.path.join(output_dir_path, 
                                     'feature-selected-datasets')
    classification_dir_path = os.path.join(output_dir_path, 
                                           'classification-results')

    os.makedirs(binary_dir_path, exist_ok=True)
    os.makedirs(oversampled_dir_path, exist_ok=True)
    os.makedirs(important_features_dir_path, exist_ok=True)
    os.makedirs(feature_selected_dir_path, exist_ok=True)
    os.makedirs(classification_dir_path, exist_ok=True)

    paths = (
        binary_dir_path, 
        oversampled_dir_path, 
        important_features_dir_path,
        feature_selected_dir_path, 
        classification_dir_path
    )
    return paths


def create_binary_datasets(multilabel_dir_path, binary_dir_path):
    '''
    Uses the create_binary_datasets.py script to convert
    a directory of mutltilabel datasets into binary
    datasets.

    Parameters
    ----------
    multilabel_dir_path: string
        Path to the directory containing multilabel datasets.
    binary_dir_path: string
        Path to save binary datasets.

    Returns
    -------
    None
    '''

    mot.create_binary_datasets.main(['-i', multilabel_dir_path, 
                                     '-o', binary_dir_path])

def create_oversampled_datasets(binary_dir_path, oversampled_dir_path, 
                                multilabel_file_name):
    '''
    Oversamples the minority class in each dataset and 
    creates train and test datasets.

    Parameters
    ----------
    binary_dir_path: string
        Path to the directory containing the binary
        dataset files.
    oversampled_dir_path: string
        Path to the directory where the resulting
        train and test files will be stored.
    multilabel_file_name: string
        The name of the multilabel dataset that
        the selected binary datasets were derived from.

    Returns
    -------
    None
    '''

    for binary_file_name in os.listdir(binary_dir_path):
        if multilabel_file_name in binary_file_name:
            binary_file_path = os.path.join(binary_dir_path, binary_file_name)
            mot.smote.main(['-i', binary_file_path, 
                             '-o', oversampled_dir_path])

def rank_features(oversampled_dir_path, important_features_dir_path,
                  multilabel_file_name, weka_path, classes_dir_path, 
                  gain_ratio_path):
    '''
    Ranks features according to their Information Gain 
    Ratio, and saves the ordered list as a csv file.

    Parameters
    ----------
    oversampled_dir_path: string
        The path to the directory containing the 
        oversampled datasets.
    important_features_dir_path: string 
        The path to the directory where the 
        important features will be saved.
    multilabel_file_name: string
        The name of the multilabel dataset that
        the selected binary datasets were derived from.

    Returns
    -------
    None
    '''

    gain_ratio_classpath = '.:' + weka_path + ':' + classes_dir_path
    for file_name in tqdm(os.listdir(oversampled_dir_path)):
        if multilabel_file_name in file_name and 'train' in file_name:
            dataset_path = os.path.join(oversampled_dir_path, file_name)
            subprocess.run([
                'javac',
                '-cp',
                gain_ratio_classpath,
                gain_ratio_path,
                '-d',
                classes_dir_path
            ])

            subprocess.run([
                'java',
                '-cp',
                gain_ratio_classpath,
                'GainRatio',
                dataset_path,
                important_features_dir_path
            ])

def feature_selection(important_features_dir_path, oversampled_dir_path,
                      feature_selected_dir_path):
    '''
    Creates new training and testing datasets using the
    top ranked features.

    Parameters
    ----------

    Returns
    -------
    None
    '''

    for file_name in tqdm(os.listdir(important_features_dir_path)):
        features_path = os.path.join(important_features_dir_path, file_name)
        tissue_loci = re.sub('_features_gain_ratio', '', file_name)
        train_path = os.path.join(oversampled_dir_path, tissue_loci)
        test_path = re.sub('train', 'test', train_path)

        mot.gain_ratio_feature_selection.main([
            '--features', features_path,
            '--train', train_path,
            '--test', test_path,
            '-o', feature_selected_dir_path
        ])

def random_forest(feature_selected_dir_path, classification_dir_path):
    '''
    Classifies the feature selected datasets using a Random
    Forest classifier and saves the results.

    Parameters
    ----------
    feature_selected_dir_path: string
        Path to the directory where the 
        feature selected datasets are stored.
    classification_dir_path: string
        Path to the directory where the
        classification results will be saved.

    Returns
    -------
    None
    '''

    for file_name in os.listdir(feature_selected_dir_path):
        if 'train' in file_name:
            test_file_name = re.sub('train', 'test', file_name)
            mot.rf_binary.main([
                '--train',
                os.path.join(feature_selected_dir_path, file_name),
                '--test',
                os.path.join(feature_selected_dir_path, test_file_name),
                '-o',
                classification_dir_path
            ])

def main(cli_args=sys.argv[1:]):
    args = parse_cli(cli_args)
    input_dir_path = args['input']
    output_dir_path = args['output']
    weka_path = args['weka_jar']
    classes_dir_path = args['classes_dir']
    gain_ratio_path = args['java_src']

    paths = create_directories(output_dir_path)
    binary_dir_path = paths[0]
    oversampled_dir_path = paths[1]
    important_features_dir_path = paths[2]
    feature_selected_dir_path = paths[3]
    classification_dir_path = paths[4]

    print('Creating Binary Datasets ...')
    create_binary_datasets(input_dir_path, binary_dir_path)
    for file_name in os.listdir(input_dir_path):
        multilabel_file_name = file_name[:-4]
        print(multilabel_file_name)
        create_oversampled_datasets(binary_dir_path, oversampled_dir_path, 
                                    multilabel_file_name)
        rank_features(oversampled_dir_path, important_features_dir_path,
                      multilabel_file_name, weka_path, classes_dir_path, 
                      gain_ratio_path)

    print('Feature Selection')
    feature_selection(important_features_dir_path, oversampled_dir_path,
                      feature_selected_dir_path)
    print('Random Forest')
    random_forest(feature_selected_dir_path, classification_dir_path)

if __name__ == '__main__':
    main()
