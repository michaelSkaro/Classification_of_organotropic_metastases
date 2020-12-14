# Classification of organotropic metastases

## Metastatic Tropism Overview

This repository is the code base for the classification of organotropic metastases. Transcriptomic profiles of 7,011 cancer patients in the TCGA database were used to classify and analyze the seeding location of primary tumors. The sequencing data and all clinicopatholigic reports for all of these patients were publicly available for bulk data mining through TCGA Biolinks.

## Dependencies
### Python
- **Version:** Python >= 3.5
- **Packages:** All required Python packages are listed in the [requirements.txt](./requirements.txt).
### Java
- **Version:** Java >= 8
- **Packages:**
  - [Weka](https://www.cs.waikato.ac.nz/ml/weka/index.html) >= 3.8.3
### R
- **Version:** >= 4.0
- **Packages:** All session.info() R packages are listed at the bottom of [Gene_Set_enrichment_and_semantic_analysis.R](feature-recapture/Gene_Set_enrichment_and_semantic_analysis.R).

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
