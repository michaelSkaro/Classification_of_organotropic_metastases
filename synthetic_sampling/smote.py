'''
SMOTE implementation for oversampling a minority class
for binary classification.
'''

import os
import sys
import csv
import argparse

import numpy as np
from sklearn.neighbors import NearestNeighbors


def parse_cli(args):
    '''
    Parses cli arguments.

    Parameters
    ----------
    args: list
        List of argument-value pairs.

    Returns
    -------
    parsed_args: dict
        Dict of parsed argument-value pairs.
    '''

    parser = argparse.ArgumentParser(
        description='SMOTE implementation for oversampling a minority class.'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Path to dataset (.csv)')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to the directory where the resulting' \
                              + ' datasets will be saved.')

    return vars(parser.parse_args(args))


def smote(minority_instances, sample_size, n_neighbors=5):
    '''
    Implementation of the SMOTE oversampling algorithm.

    Parameters
    ----------
    minority_instances: NumPy Array
        The instances from the dataset that belong the
        minority class.
    sample_size: int
        The number of synthetic samples to generate.
    n_neighbors: int
        The number of neighbors used during the
        oversampling process.

    Returns
    -------
    synthetic_samples: NumPy Array(sample_size, n_features)
        The synthetic samples of the minority class.
    '''

    neigh = NearestNeighbors(n_neighbors=5)
    neigh.fit(minority_instances)
    neigh_distances, neigh_ind = neigh.kneighbors(minority_instances)
    synthetic_samples = []
    for i in range(sample_size):
        sample_idx = np.random.randint(len(minority_instances))
        neigh_idx = neigh_ind[sample_idx][np.random.randint(n_neighbors)]
        diff = \
            minority_instances[neigh_idx] - minority_instances[sample_idx]
        weight = np.random.rand()
        while weight == 0:
            weight = np.random.rand()   

        weighted_diff = diff * weight
        synthetic_sample = minority_instances[sample_idx] + weighted_diff
        synthetic_samples.append(synthetic_sample)

    synthetic_samples = np.array(synthetic_samples, dtype=np.float)
    return synthetic_samples


def create_train_test_datasets(synthetic_samples, x_majority, x_minority,
                                  attr_names, majority_label, minority_label,
                                  real_minority_split=0.8):
    '''
    Creates train and test datasets using real and
    synthetic data.

    Parameters
    ----------
    synthetic_samples: NumPy Array(sample_size, n_features)
        The synthetic samples of the minority class.
    x_majority: NumPy Array
        The instances of the dataset that belong the
        majority class.
    x_minority: NumPy Array
        The instances of the dataset that belong the
        minority class.
    attr_names: list of strings
        Names of the attributes of the data.
    majority_label: string
        Name of the label describing the majority class.
    minority_label: string
        Name of the label describing the minority class.
    real_minority_split: float
        Percentage of the minority class to include in
        the training dataset.

    Returns
    -------
    train_dataset, test_dataset: NumPy Arrays
        Training dataset that includes synthetic samples,
        and a testing dataset that only includes samples
        from the real dataset.
    '''

    synthetic_samples = np.column_stack([
        synthetic_samples, 
        [minority_label for i in range(synthetic_samples.shape[0])]
    ])

    x_majority = np.column_stack([
        x_majority,
        [majority_label for i in range(x_majority.shape[0])]
    ])

    x_minority = np.column_stack([
        x_minority,
        [minority_label for i in range(x_minority.shape[0])]
    ])
    
    majority_inds = np.arange(x_majority.shape[0])
    np.random.shuffle(majority_inds)
    minority_inds = np.arange(x_minority.shape[0])
    np.random.shuffle(minority_inds)

    train_dataset = np.vstack([
        attr_names[2:], 
        synthetic_samples,
        x_majority[majority_inds[:synthetic_samples.shape[0]]],
        
    ])

    test_dataset = np.vstack([
        attr_names[2:],
        x_minority,
        x_majority[synthetic_samples.shape[0]:]
    ])

    return train_dataset, test_dataset


def write_datasets(train_dataset, test_dataset, file_name, outdir_path):
    '''
    Writes the training and testing datasets to disk.

    Parameters
    ----------
    train_dataset, test_dataset: NumPy Arrays
        Training dataset that includes synthetic samples,
        and a testing dataset that only includes samples
        from the real dataset.
    file_name: string
        Name of the input dataset.
    outdir_path: string
        Path to the directory to save the datasets.

    Returns
    -------
    None
    '''

    train_outfile_path = os.path.join(outdir_path, file_name + '_train.csv')
    test_outfile_path = os.path.join(outdir_path, file_name + '_test.csv')

    with open(train_outfile_path, 'w+') as fp:
        writer = csv.writer(fp)
        for row in train_dataset:
            writer.writerow(row)

    with open(test_outfile_path, 'w+') as fp:
        writer = csv.writer(fp)
        for row in test_dataset:
            writer.writerow(row)


def main():
    args = parse_cli(sys.argv[1:])
    with open(args['input']) as fp:
        reader = csv.reader(fp)
        data = list(reader)

    data = np.array(data)
    attr_names = data[0]
    x = data[1:, 2:-1].astype(np.float)
    y = data[1:, -1]
    classes = np.unique(y)

    class_counts = np.empty(len(classes)).astype(np.int)
    for i in range(len(class_counts)):
        class_counts[i] = np.sum(y == classes[i])

    majority_idx = np.where(np.max(class_counts) == class_counts)[0][0]
    if majority_idx == 0:
        minority_idx = 1
    else:
        minority_idx = 0

    x_majority = []
    x_minority = []
    for i, class_name in enumerate(y):
        if class_name == classes[majority_idx]:
            x_majority.append(x[i])
        else:
            x_minority.append(x[i])

    x_majority = np.array(x_majority, dtype=np.float)
    x_minority = np.array(x_minority, dtype=np.float)
    #print('Majority Instances:', x_majority.shape, 
    #      'Minority Instances:', x_minority.shape)
    sample_size = int((len(x_majority) - len(x_minority)) * .8)
    synthetic_samples = smote(x_minority, sample_size)

    train_dataset, test_dataset = create_train_test_datasets(
        synthetic_samples, 
        x_majority, 
        x_minority,
        attr_names, 
        classes[majority_idx], 
        classes[minority_idx]
    )

    file_name = os.path.split(args['input'])[1].split('.')[0]
    write_datasets(train_dataset, test_dataset, file_name, args['output'])

 
if __name__ == '__main__':
    main()
