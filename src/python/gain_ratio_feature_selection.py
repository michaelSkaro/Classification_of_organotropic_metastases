'''
Creates datasets using the features with the highest gain
ratio scores.
'''

import os
import re
import csv
import sys
import argparse

import numpy as np

def parse_cli(args):
    '''
    Parses command line arguments.
    '''

    parser = argparse.ArgumentParser(
        description='Creates datasets using the features with the highest' 
                     + ' gain ratio scores.'
    )
    parser.add_argument('--features', required=True, 
                        help='Path to features (.csv).')
    parser.add_argument('--train', required=True, 
                        help='Path to training dataset (.csv).')
    parser.add_argument('--test', required=True, 
                        help='Path to testing dataset (.csv).')
    parser.add_argument('-o', '--output', required=True, 
                        help='Path to output directory.')
    return vars(parser.parse_args(args))

def main():

    args = parse_cli(sys.argv[1:])
    with open(args['features']) as fp:
        reader = csv.reader(fp)
        features = list(reader)

    with open(args['train']) as fp:
        reader = csv.reader(fp)
        train_data = list(reader)

    with open(args['test']) as fp:
        reader = csv.reader(fp)
        test_data = list(reader)

    features = np.array(features)[1:, 0]
    train_data = np.array(train_data)
    test_data = np.array(test_data)

    reduced_train_data = train_data[:,-1]
    reduced_test_data = test_data[:,-1]
    for feature in features:
        idx = np.where(train_data[0] == feature)[0][0]
        reduced_train_data = np.column_stack([
            train_data[:, idx], 
            reduced_train_data
        ])
        reduced_test_data = np.column_stack([
            test_data[:, idx], 
            reduced_test_data
        ])

    
    train_file_name = os.path.split(args['train'])[1]
    output_train_name = re.sub('_train', '_feature_selected_train', train_file_name)
    output_test_name = re.sub('train', 'test', output_train_name)
    output_train_path = os.path.join(args['output'], output_train_name)
    output_test_path = os.path.join(args['output'], output_test_name)
    with open(output_train_path, 'w+') as fp:
        writer = csv.writer(fp)
        for row in reduced_train_data:
            writer.writerow(row)

    with open(output_test_path, 'w+') as fp:
        writer = csv.writer(fp)
        for row in reduced_test_data:
            writer.writerow(row)
        

if __name__ == '__main__':
    main()
