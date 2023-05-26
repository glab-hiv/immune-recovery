# immune-recovery
## CD4 T cell receptor hierarchies are stable and independent of HIV-mediated dysregulation of immune homeostasis
***Link publication here***

## CyTOF Data
A 31-marker mass cytometry panel was used to examine all major PBMC lineages and specifically CD4 and CD8 T cell memory dynamics in people with HIV (PWH) who are durably ART suppressed for an average of 6.7 years (LT-ART, n=10) and PWH in the first 500 days following ART initiation (A5248, n=10). The panel also includes markers of activation (HLA-DR, CD38, CCR5), activation/exhaustion (PD-1), proliferation (Ki67), survival (Bcl-2) and long-lived memory (CD127).

You can download the publicly available preprocessed CyTOF files (FCS format) or preprocessed downsampled data from the unsupervised analysis (h5ad format) from the [Zenodo](https://zenodo.org/record/7495836) repository.

### Cohort Summary 
<p align="center">
  <img src="overview.png" width="75%" height ="75%" />
</p>

## Installation
You can clone the git repository by, 
```
git clone https://github.com/glab-hiv/immune-recovery.git
```

## Elastic Net Regression
### Dependencies
* Python >= 3.6, anndata 0.7.6, scanpy 1.8.1, numpy 1.19.5, scipy 1.7.1, leidenalg 0.8.7, flowkit 0.9.1, umap-learn 0.5.2, scikit-learn 0.24.1, phenograph 1.5.7

### Author
Jolene Ranek, <ranekj@live.unc.edu>

### Description
To ascertain immune cell phenotype and functional changes associated with initiation of antiretroviral therapy or durable suppression, a predictive elastic net regression model was trained on immune features that described the cell type composition and cell type signaling activity of profiled samples over time. Python scripts for preprocessing (preprocess.py, subsample.py), clustering (cluster.ipynb), and elastic net regression (run_elasticnet.ipynb) are described below. 

* `preprocess.py` - Creates a preprocessed __.h5ad__ data object for each patient cohort. This script parses FCS files from the Zenodo data folder and performs preprocessing by arcsinh transforming data with a cofactor of 5. 
* `subsample.py` - Downsamples data (2500 cells per patient sample) using [Kernel Herding sketching](https://dl.acm.org/doi/abs/10.1145/3535508.3545539).
* `cluster.ipynb` - Performs unsupervised meta-clustering on the downsampled data using the [PhenoGraph algorithm](https://pubmed.ncbi.nlm.nih.gov/26095251/).
* `run_elasticnet.ipynb` - Engineers descriptive mass cytometry immune features from unsupervised clusters (e.g. cell population frequency) and performs [elastic net regression](https://www.jstor.org/stable/3647580) to identify changes in unsupervised cluster-derived features over time. 

## Quasi-Binomial and Gamma Fixed-Effects Regression
### Dependencies
* R >= 4.2.2, ggplot2>=3.4.0, dplyr>=1.0.10, margins 0.3.26, writexl 1.4.1

### Author
Ann Marie Weideman, <anndo1@umbc.edu> (preferred) or <anndo1@live.unc.edu>

### Description
To examine longitudinal changes in the proportion of different cell populations (e.g., CD4 cells) expressing certain biomarkers (e.g., CCR5+ or a combination of biomarkers CCR7+CD45RA+), [quasi-binomial fixed-effects regression models](https://books.google.com/books/about/Generalized_Linear_Models_Second_Edition.html?id=h9kFH2_FfBkC) were fit separately for data from each cell population and each cohort (LT-ART and A5248) using the data in the Zenodo data folder. 

Each R script described below plots the data to check model assumptions (e.g., linearity to the logit), fits a separate model for each cell population and each cohort (LT-ART and A5248), and outputs the slope estimates and corresponding confidence assessments with accompanying interpretation. 

* `QuasiBinom_FEM_CD4.r` quasi-binomial fixed-effects regression model for CD4 T cell memory dynamics. 
* `QuasiBinom_FEM_CD8.r` quasi-binomial fixed-effects regression model for CD8 T cell memory dynamics.

Longitudinal changes in proteins expressed on the cell surface and intracellularly as measured by mean signal intensity (MSI) were also examined. Since values of MSI are positive and continuous, a gamma fixed-effects regression with log link was fit separately to the data from each expressed protein and each cohort (LT-ART and A5248) using the data in the Zenodo data folder, and deviance residuals were inspected for normality. 

Each R script described below fits a separate model for each cell population and each cohort (LT-ART and A5248), checks the deviance residuals for normality, and outputs the slope estimates and corresponding confidence assessments with accompanying interpretation.

* `Gamma_FEM_CD4.r` gamma fixed-effects regression model for CD4 T cell memory dynamics. 
* `Gamma_FEM_CD8.r` gamma fixed-effects regression model for CD8 T cell memory dynamics.
