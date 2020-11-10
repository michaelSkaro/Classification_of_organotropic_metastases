'''
Creates new dataset of equal cancer and normal class
proportions from TCGA data.
'''

import os
import re
import csv
import argparse

import joblib
import numpy as np
from tqdm import tqdm

def create_new_dataset(input_path, label, output_path):
    '''
    Construct new dataset based on the input file.

    Parameters
    ----------
    input_path: string
        Path to the TCGA dataset (.csv).
    label: string
        Label of classes with the smallest number of 
        instances in the binary dataset.
    output_path: string
        Path to an output directory to save the newly 
        created dataset.

    Returns
    ----------
    NoneType object
    '''
    
    print(input_path)
    with open(input_path, 'r') as fp:
        reader = csv.reader(fp, delimiter=',')
        data = list(reader)
        data = np.array(data)

    data[0,0] = " "
    smaller_class = [x for x in data if x[-1] == label]
    larger_class = [x for x in data if x[-1] != label]
    np.random.shuffle(larger_class)
    
    file_name = os.path.split(input_path)[-1]
    out_file = os.path.join(
        output_path, 
        re.sub('learning', 'balanced', file_name)
    )
    with open(out_file, 'w') as fp:
        writer = csv.writer(fp, delimiter=',')
        writer.writerow(data[0])
        for i in range(len(smaller_class)):
            writer.writerow(larger_class[i])
            writer.writerow(smaller_class[i])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates new dataset of \
                                                  equal cancer and normal \
                                                  class proportions from \
                                                  binary TCGA data.')
    parser.add_argument('-i', '--input', required=True, 
                        help='Input dataset path or a directory \
                              of datasets.')
    parser.add_argument('-l', '--label', required=True,
                        help='Label of the class with the least instances.')
    parser.add_argument('-o', '--output', default=os.getcwd(), 
                        help='Path to an output directory to save the dataset.')
    parser.add_argument('-n', '--n_jobs', default=-1, type=int,
                        help='Max number of concurrent jobs.')
    args = vars(parser.parse_args())

    if os.path.isdir(args['input']):
        inputs = []
        for file_name in os.listdir(args['input']):
            if 'learning' in file_name:
                inputs.append(os.path.join(args['input'], file_name))
    else:
        inputs = [args['input']]

    joblib.Parallel(n_jobs=args['n_jobs'])(
        joblib.delayed(create_new_dataset)(
            input_file, args['label'], args['output']
        )
        for input_file in tqdm(inputs)
    )
