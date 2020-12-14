'''
Creates binary datasets from multilabel data.
'''

import os
import csv
import sys
import argparse

import numpy as np

def parse_cli(args):
    '''
    '''

    parser = argparse.ArgumentParser(
        description='Creates binary datasets from multilabel data.'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Path to directory of input datasets.')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to directory to save datasets.')

    return vars(parser.parse_args(args))

def main(cli_args=sys.argv[1:]):

    args = parse_cli(cli_args)
    indir_path = args['input']
    outdir_path = args['output']
    for file_name in os.listdir(indir_path): #tqdm
        data_path = os.path.join(indir_path, file_name)
        with open(data_path) as fp:
            reader = csv.reader(fp)
            data = list(reader)
                
        data =  np.array(data)
        label_indices = [-1 if 'ENSG' in x else i for i, x in enumerate(data[0])] 
        labels_starting_index = np.min([x for x in label_indices if x > 1])
        num_of_labels = len(data[0]) - labels_starting_index

        for i in range(num_of_labels):
            label_index = labels_starting_index + i
            outfile_name = file_name[:-4] + '_' + data[0, label_index] + '.csv'
            outfile_path = os.path.join(outdir_path, outfile_name)
            with open(outfile_path, 'w+') as fp:
                writer = csv.writer(fp, delimiter=',')
                for j, row in enumerate(data):
                    output = list(row[:labels_starting_index])
                    if j == 0:
                        output.append(row[label_index])
                    else:
                        if row[label_index] == '1':
                            val = 'Positive'
                        else:
                            val = 'Negative'
                        output.append(val)
                    
                    writer.writerow(output)

if __name__ == '__main__':
    main()

