# Classification of organotropic metastases

## Metastatic Tropism Overview

This repository is the code base for the classification of organotropic metastases. Transcriptomic profiles of 7,011 cancer patients in the TCGA database were used to classify and analyze the seeding location of primary tumors. The sequencing data and all clinicopatholigic reports for all of these patients were publicly available for bulk data mining through TCGA Biolinks.

## Installation
We utilized multiple programming languages (i.e. Java, Python, and R) to construct learning models and to perform biological analyses. As a result, this created many dependencie so we have provided two different ways to install our framework -- a manual installation and a docker installation.

**Docker Installation**

The docker image for this project can be pulled from the online Docker Hub [repository](https://hub.docker.com/r/marcdh3/mot) or can be built using the Dockerfile included in the base directory of this project.

To pull the image from the Docker Hub repo, run the following command:
```
docker pull mskaro1/mot
```

To build the image using the Dockerfile, run the following command in the base directory of this project:
```
docker build --tag mskaro1/mot .
```

**Manual Installation**

For those seeking to manually install the project, all of the following dependencies must be satisfied prior to attempting the installation:

Python
- Version: Python >= 3.5
- Packages: All required Python packages are listed in the [requirements.txt](./requirements.txt).

Java
- Version: Java >= 8
- Packages:
  - [Weka](https://www.cs.waikato.ac.nz/ml/weka/index.html) >= 3.8.3 (Note: The jar file is already included [here](./lib/weka.jar))
  
R
- Version: >= 4.0
- Packages: All session.info() R packages are listed at the bottom of [Gene_Set_enrichment_and_semantic_analysis.R](feature-recapture/Gene_Set_enrichment_and_semantic_analysis.R).

Once all of the required depenedencies are satisfied, run the following command in the base directory to install the project as a python package:
```
pip install .
```

## Metastatic Classification Demo:
We have provided a sample dataset of TCGA data to demonstrate the effectiveness of our metastatic classification approach. Our sample [dataset](./samples/metastasis-demo/TCGA-COAD_metastatic_data_RNAseq.csv) is of Colon Adenocarcinoma (COAD) tumors that metastasiszed to the colon, liver, or lung.

**Recommneded: Docker Approach**
```
docker run --rm -it -v <output-directory>:/demo-outputs mskaro1/mot
```
Note: `<output-directory>` should be replaced with the path of a directory on the user's local machine, and it is where the outputs of the demo will be stored.
  
**Manual Approach**
```
 python3 -m mot.metastasis_pipeline -i ./samples/metastasis-demo/ -o <output-directory> -w ./lib/weka.jar -c ./classes -j /src/GainRatio.java
```
Note: This command should be run in the base directory of the project, and `<output-directory>` should be replaced with the path to a directory for the outputs to be stored.
  
**Demo Outputs**
```
├── <output-directory>/
│   ├── binary-datasets/
│   ├── oversampled-datasets/
│   ├── important-features/
|   ├── feature-selected-datasets/
|   ├── classification-results/
```
- binary-datasets: The multilabel COAD dataset is split into multiple binary datasets, and the binary datasets are stored in this directory. 
- oversampled-datasets: The training and testing data generated from the binary datasets. The training data uses synthetic data generated by the SMOTE algorithm, while the testing data uses only real TCGA data.
- important-features: The top 1000 features (i.e. genes) of each training dataset ranked by their information gain ratio score.
- feature-selected-datasets: The training and testing datasets that only contain the top 1000 selected features.
- classification-results: Directory contains the classification results of our Random Forest model on the feature-selected datasets.

## General Usage:
The entire metastatic pipeline can be ran using the [metastasis_pipeline script](./src/metastasis_pipeline.py). This script is callable from the command line interface using the following command: 

```python -m mot.metastasis_pipeline``` 

The `-h` flag to understand all available options. 

Additionally, each component of the pipeline can be called individually from the command line. For more information read our [wiki](./) for a breakdown of each script's role in the pipeleine.

**Note:** For those seeking to use the docker image to interact with our framework, run the following command to gain access to the shell of the docker image:
```
docker run --rm -it --entrypoint="" mskaro1/mot bash
```

## Classification Stages
1. [Cancer vs. Normal Classification](cancer-vs-normal-classification)
    - [undersample.py](cancer-vs-normal-classification/undersample.py): Creates a new dataset of equal cancer and normal class proportions from TCGA data.
    - [RandomForestClassifier.java](cancer-vs-normal-classification/RandomForestClassifier.java): Constructs a Random Forest, performs 10-Fold cross-validation, and saves the results.
2. [Cancer vs. Cancer Classification](cancer-vs-cancer-classification)
    - [multiclass_cancer_sampling.py](cancer-vs-cancer-classification/multiclass_cancer_sampling.py): Generates either a single balanced, multiclass dataset containing samples from each cancer tissue type or multiple balanced, binary datasets of all pairwise cancer type combinations.
    - [RandomForestClassifier.java](cancer-vs-normal-classification/RandomForestClassifier.java): Constructs a Random Forest, performs 10-Fold cross-validation, and saves the results. (The same code used during the 1st Stage was used during this stage as well.)
3. [Metastases Classification](metastases-classification)
    - [create_binary_datasets.py](metastases-classification/create_binary_datasets.py): Creates binary datasets from multilabel data.
    - [smote.py](synthetic-sampling/smote.py): SMOTE implementation for oversampling a minority class for binary classification.
    - [GainRatio.java](metastases-classification/GainRatio.java): Ranks the attributes of a dataset by their gain ratio.
    - [gain_ratio_feature_selection.py](metastases-classification/gain_ratio_feature_selection.py): Creates datasets using the features with the highest gain ratio scores.
    - [feature_sapce_plot.py](metastases-classification/feature_sapce_plot.py): Multidimensional Scaling plots of feature space prior too and following feature selection
    - [rf_binary.py](metastases-classification/rf_binary.py): Random Forest implementation with code for hyper-parameter tuning.
    - [Viz_MOT.ipynb](metastases-classification/Visualization_for_MOT.ipynb): Classification of positive and negative progression to seeding locations in 6 cancer types
    - [metastasis_pipeline.py](metastases-classification/metastasis_pipeline.py): A supervisory script to run all data through the metastatic loci classificaiton pipeline.
4.  [Feature recapture and analysis](feature-recapture)
    - [Enriched_features_Fisher_recap.R](feature-recapture/Enriched_features_Fisher_recap.R): Statistical analysis of recaptured features in independent lists.
    - [Gene_Set_enrichment_and_semantic_analysis.R](feature-recapture/Gene_Set_enrichment_and_semantic_analysis.R): Gene set enrichment analysis in GO Database. Sematic analysis of GO biological processes.  
    
## Reviewing:

We are making this pipeline into a Docker image. Stay tuned!
Thanks!
