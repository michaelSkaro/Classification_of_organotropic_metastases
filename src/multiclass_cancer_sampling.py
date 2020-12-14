'''
Generates either a single balanced, multiclass dataset
containing samples from each cancer tissue type, or 
multiple balanced, binary datasets of all pairwise
cancer type combinations.
'''

import os
import csv
import sys
import argparse

import numpy as np
from tqdm import tqdm

def all_cancers(tcga_dir, outfile, count, label, use_label=False):
    '''
    Creates a dataset from the TCGA files.

    Parameters
    ----------
    tcga_dir: string
        Path to the directory containing TCGA read counts.
    outfile: string
        Path to save the dataset.
    count: int
        The number of instances of each class type to be in
        the dataset.
    label: string
        Specific class in each dataset to sample instances 
        from.
    use_label: boolean

    Returns
    ----------
    NoneType object
    '''
    
    tcga_files = []
    for file_name in os.listdir(tcga_dir):
        if 'learning' in file_name:
            tcga_files.append(os.path.join(tcga_dir, file_name))
    
    with open(outfile, 'a+') as fp: #w
        writer = csv.writer(fp, delimiter=',')
        for i, file_path in enumerate(tqdm(tcga_files)):
            with open(file_path, 'r') as fp:
                reader = csv.reader(fp, delimiter=',')
                data = list(reader)
                data = np.array(data)
                data[0,0] = ' '

            cancers = data[data[:,-1] == label]
            np.random.shuffle(cancers)
            file_name = os.path.split(file_path)[-1]

            if not use_label:
                cancers[:,-1] = file_name.split('.')[0].split('_')[0]

            if i == 0:
                writer.writerow(data[0])

            for instance in cancers[:count]:
                writer.writerow(instance)
    
def cancer_pairs(tcga_dir, out_dir, count, label):
    '''
    Generates all datasets of pairwise combiniations of 
    cancer types.

    Parameters
    ----------
    tcga_dir: string
        Path to the directory containing TCGA read counts.
    out_dir: string
        Path to directory to save datasets.
    count: int
        The number of instances of each class type to be in
        the dataset.
    label: string
        Specific class in each dataset to sample instances 
        from.

    Returns
    ----------
    NoneType object
    '''

    #Determine file combinations (n choose 2)
    files = [x for x in os.listdir(tcga_dir) if 'learning' in x]
    pairs = []
    for x in files:
        for y in files:
            if (y, x) not in pairs and x != y:
                pairs.append((x, y))

    for row in tqdm(pairs):
        file_name = row[0].split('_')[0] + '_' + row[1].split('_')[0] + '.csv'
        outfile = os.path.join(out_dir, file_name)
        with open(outfile, 'w') as fp:
            writer = csv.writer(fp, delimiter=',')
            for i, input_file in enumerate(row):
                input_path = os.path.join(tcga_dir, input_file)
                with open(input_path, 'r') as file_reference:
                    reader = csv.reader(file_reference, delimiter=',')
                    data = list(reader)
                    data = np.array(data)
                    data[0,0] = ' '
                    cancers = data[data[:,-1] == label]
                    cancers[:,-1] = input_file.split('.')[0].split('_')[0]
                    np.random.shuffle(cancers)
                    if i == 0:
                        writer.writerow(data[0])
                    for instance in cancers[:count]:
                        writer.writerow(instance)

def comprehensive(input_dir, outfile, counts):
    '''
    Extracts the labels necessary to distinguish 
    cancer tissue samples from normal
    tissue samples (i.e. Primary_Tumor and 
    Solid_Tissue_Normal). Then the all_cancers function
    is called using the labels to sample each dataset
    for the appropriate instances.

    Parameters
    ----------
    input_dir: string
        Path to the directory containing TCGA read counts.
    outfile: string
        Path to save the dataset.
    counts: int
        The number of instances of each class type to be in
        the dataset.

    Returns
    -------
    None
    '''

    #Determine the labels
    filename = os.listdir(input_dir)[0]
    with open(os.path.join(input_dir, filename)) as fp:
        reader = csv.reader(fp)
        data = np.array(list(reader))

    labels = np.unique(data[1:,-1])
    for label in labels:
        print(label)
        all_cancers(input_dir, outfile, counts, label, True)
    


def parse_cli(args):
    '''
    Parse command line arguments.

    Parameters
    ----------
    args: list of strings
        argument options and their values

    Returns
    ----------
    args: dict
        A dict that argument values can be referenced in 
        by supplying the option or flag name as the key.
    '''

    parser = argparse.ArgumentParser(description='Generate balanced datasets'
                                                 + ' of different cancer'
                                                 + ' types.')
    parser.add_argument('-i', '--input', required=True, 
                        help='Input directory of TCGA files (.csv)')
    parser.add_argument('-o', '--output', required=True,  
                        help='If pairs flag is used then output must be'
                              + ' a directory path. Otherwise, output must be'
                              + ' a csv file path.')
    parser.add_argument('-c', '--count', required=True, type=int, 
                        help='Number of instances per class.')
    parser.add_argument('-l', '--label', default='Primary_Tumor', 
                        help='Specific class in each dataset to sample'
                             + '  instances from.')
    parser.add_argument('--comprehensive', action='store_true', default=False,
                        help='Use this flag to generate a comprehensive'
                              + ' dataset of all tissue types and their'
                              + ' classes.')
    parser.add_argument('--pairs', action='store_true', default=False,
                        help='Use this flag to generate all pairwise'
                              + ' combinations of cancer types.')

    return vars(parser.parse_args(args))
    
def main():
    args = parse_cli(sys.argv[1:])
    if not os.path.isdir(args['input']):
        sys.exit('Input must be a directory. Use -h or --help for more'
                  + ' information.')

    if args['comprehensive']:
        comprehensive(args['input'], args['output'], args['count'])
    elif args['pairs'] and not os.path.isfile(args['output']):
        os.makedirs(args['output'], exist_ok=True)
        cancer_pairs(args['input'], args['output'], args['count'], args['label'])
    elif not args['pairs'] and not os.path.isdir(args['output']):
        all_cancers(args['input'], args['output'], args['count'], args['label'])
    else:
        sys.exit('Invalid output type. Use -h or --help for more information.')

if __name__ == '__main__':
    main()      
