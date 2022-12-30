# immune-recovery
## Placeholder for Title/Link to Publication
***
## CyTOF Data (**sent to Alexis to edit**)
Includes preprocessed CyTOF files (FCS format) from a [***]-marker mass cytometry panel to examine CD4 and CD8 T cell memory dynamics in people with HIV (PWH) who are durably (~5yrs) ART suppressed (PWH-ART, n=10) and PWH in the 1.5 years following ART initiation (n=10, ACTG5248). The panel included markers of activation (HLA-DR, CD38, CCR5), activation/exhaustion (PD-1), proliferation (Ki67), survival (Bcl-2) and long-lived memory (CD127). 

### Cohort Summary 
1. People with HIV durably suppressed on ART (PWH-ART)
   - n=10 participants
   - 3 different timepoints for each
2. People with HIV in the 1.5 years following ART initiation (ACTG5248)
   - n=10 participants
   - 12 timepoints for each at days: 0, 2, 7, 10, 14, 21, 28, 56, 84, 140, 252, 504

## Elastic Net Regression
***
## Quasi-Binomial Fixed-Effects Regression
To examine whether there was a change in leukocyte populations and subpopulations over time, a quasi-binomial fixed-effects regression with indicators for participant ID (PID) was fit separately to the data from each cohort (PWH-ART and ACTG5248) using the data in the CyTOF Data folder. The R code for these models is separated into two scripts for CD4 T cell memory dynamics (FEreg_CD4.r) and for CD8 T cell memory dynamics (FEreg_CD8.r). Each script plots the data to check model assumptions (e.g., linearity to the logit), fits a separate model for each population, and outputs the slope estimates and corresponding confidence assessments with accompanying interpretation. 
