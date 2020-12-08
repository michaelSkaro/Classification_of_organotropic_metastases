'''
Creates MDS plots of the metastatic data pre-, post-feature-selection.
'''

import os
import re
import csv

import numpy as np
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def save_dual_plots(file_name, pre_components, post_components,
                    y_pre_data, y_post_data, output_dir_path):
    '''
    Save plots that are side-by-side of the 
    pre-, post-feature selection.

    Parameters
    ----------
    file_name: string
        Name of the feature selected dataset file.
    pre_components: NumPy array (NxM)
        PCA components of the pre-feature selection data,
        where N is the number of instances in the dataset 
        and M is the number of attributes.
    post_components: NumPy array (NxM)
        PCA components of the post-feature selection data,
        where N is the number of instances in the dataset
        and M is the number of attributes.
    y_pre_data: NumPy array (N,)
        Class labels of the data, where N is the number of
        instances in the pre-feature selection dataset.
    y_post_data: NumPy array (N,)
        Class labels of the data, where N is the number of
        instances in the post-feature selection dataset.
    output_dir_path: string
        Path to the directory where the plots will be saved.

    Returns
    -------
    None
    '''

    suptitle = re.sub('TCGA-', '', file_name)
    suptitle = re.sub('metastatic_data_RNAseq_', '', suptitle)
    suptitle = re.sub('_feature_selected_test', '', suptitle)
    suptitle = re.sub('.csv', '', suptitle)
    fig = plt.figure(figsize=(10, 10))
    plt.suptitle(suptitle)

    ax = fig.add_subplot(121)
    ax.set_title('Pre-Feature-Selection')
    ax.set_ylabel('2nd Principal Component')
    ax.set_xlabel('1st Principal Component')
    sns.kdeplot(x=pre_components[:, 0], y=pre_components[:, 1], 
                hue=y_pre_data, ax=ax)

    ax = fig.add_subplot(122)
    ax.set_title('Post-Feature-Selection')
    ax.set_ylabel('2nd Principal Component')
    ax.set_xlabel('1st Principal Component')
    sns.kdeplot(x=post_components[:, 0], y=post_components[:, 1], 
                hue=y_post_data, ax=ax)

    plt.savefig(os.path.join(output_dir_path, suptitle + '.png'), dpi=1000)
    plt.close()

def save_individual_plots(file_name, pre_components, post_components, 
                          y_pre_data, y_post_data, output_dir_path):
    '''
    Save the pre-, post-feature selection plots as 
    individual figures.

    Parameters
    ----------
    file_name: string
        Name of the feature selected dataset file.
    pre_components: NumPy array (NxM)
        PCA components of the pre-feature selection data,
        where N is the number of instances in the dataset 
        and M is the number of attributes.
    post_components: NumPy array (NxM)
        PCA components of the post-feature selection data,
        where N is the number of instances in the dataset
        and M is the number of attributes.
    y_pre_data: NumPy array (N,)
        Class labels of the data, where N is the number of
        instances in the pre-feature selection dataset.
    y_post_data: NumPy array (N,)
        Class labels of the data, where N is the number of
        instances in the post-feature selection dataset.
    output_dir_path: string
        Path to the directory where the plots will be saved.

    Returns
    -------
    None
    '''

    '''
    # Swaps instances to ensure that plot curve colors 
    # are consistent. 
    temp_comp = pre_components[0]
    pre_components[0] = pre_components[23]
    pre_components[23] = temp_comp
    temp_data = y_pre_data[0]
    y_pre_data[0] = y_pre_data[23]
    y_pre_data[23] = temp_data
    '''

    temp_comp = post_components[0]
    post_components[0] = post_components[23]
    post_components[23] = temp_comp
    temp_data = y_post_data[0]
    y_post_data[0] = y_post_data[23]
    y_post_data[23] = temp_data

    sns.set_context(context='paper', font_scale=1)
    suptitle = re.sub('TCGA-', '', file_name)
    suptitle = re.sub('metastatic_data_RNAseq_', '', suptitle)
    suptitle = re.sub('_feature_selected_test', '', suptitle)
    suptitle = re.sub('.csv', '', suptitle)
    figsize = (5, 5) #(8, 8)

    fig = plt.figure(figsize=figsize)
    plt.suptitle(suptitle, fontweight='bold', fontsize=15)
    ax = fig.add_subplot(111)
    ax.set_title('Pre-Feature-Selection', fontsize=15)
    ax.set_ylabel('2nd Principal Component', fontsize=15)
    ax.set_xlabel('1st Principal Component', fontsize=15)
    sns.kdeplot(x=pre_components[:, 0], y=pre_components[:, 1], 
                hue=y_pre_data, ax=ax)
    #plt.legend(prop={'size':13}, loc='lower right')
    plt.savefig(os.path.join(output_dir_path, suptitle + '_pre_feature_selection.png'),)# dpi=1000)
    plt.close()
    
    fig = plt.figure(figsize=figsize)
    plt.suptitle(suptitle, fontweight='bold', fontsize=15)
    ax = fig.add_subplot(111)
    ax.set_title('Post-Feature-Selection', fontsize=15)
    ax.set_ylabel('2nd Principal Component', fontsize=15)
    ax.set_xlabel('1st Principal Component', fontsize=15)
    sns.kdeplot(x=post_components[:, 0], y=post_components[:, 1], 
                hue=y_post_data, ax=ax)
    #plt.yticks(np.arange(-2e6, 6e6, 1e6))
    plt.savefig(os.path.join(output_dir_path, suptitle + '_post_feature_selection.png'),)# dpi=1000)
    plt.close()

def main():
    '''
    Pre-feature-selected dataset: oversampled
    Post-feature-selected dataset: feature-selected-dataset

    Use test files.
    '''

    sns.set_context('paper')
    pre_dir_path = '/home/marcus/Desktop/pipeline-outputs/oversampled-datasets'
    post_dir_path = '/home/marcus/Desktop/pipeline-outputs/feature-selected-datasets'
    output_dir_path = '/home/marcus/Desktop/metastatic-pca-plots'

    for file_name in tqdm(os.listdir(post_dir_path)):
        if 'test' in file_name:
            pre_file_name = re.sub('_feature_selected', '', file_name)
            
            with open(os.path.join(pre_dir_path, pre_file_name)) as fp:
                reader = csv.reader(fp)
                data = list(reader)

            data = np.array(data)
            x_pre_data = data[1:, :-1].astype(np.float)
            y_pre_data = data[1:, -1]
            pre_components = PCA().fit_transform(x_pre_data)

            with open(os.path.join(post_dir_path, file_name)) as fp:
                reader = csv.reader(fp)
                data = list(reader)

            data = np.array(data)
            x_post_data = data[1:, :-1].astype(np.float)
            y_post_data = data[1:, -1]
            post_components = PCA().fit_transform(x_post_data)

            #save_dual_plots(file_name, pre_components, post_components, 
            #                y_pre_data, y_post_data, output_dir_path)       
            save_individual_plots(file_name, pre_components, post_components, 
                                  y_pre_data, y_post_data, output_dir_path)
            
    
    return 0

if __name__ == '__main__':
    main()
