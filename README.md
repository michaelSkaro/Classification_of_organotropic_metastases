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

## Classification Stages
"Summary of all the classification stages"

- [Cancer vs. Normal Classification](cancer-vs-normal)
- [Multiclass Cancer Classification](multiclass-cancer-classification)
- [Metastases Classification](metastases-classification)
