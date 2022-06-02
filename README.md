Algorithm
-------
Replicability in Cancer Omics Data Analysis

Maintainer
-------
Jiping Wang   <jiping.wang@yale.edu>


Publication
-------
Jiping Wang, Hongmin Liang, Qingzhao Zhang, Shuangge Ma (2022). Replicability in Cancer Omics Data Analysis: Measures and Empirical Explorations. Manuscript


Usage
-------
 The following .R files are submitted.
 
(Data Processing: TCGA data pre-processing to generate the final data for survival analysis)

1.	1-Data Clean and Query Code-FPKM/TPM.R: This file includes steps processing the TCGA raw FPKM data downloaded directly from the NCI Genomic Data Commons Data Portal (https://portal.gdc.cancer.gov/) into TPM format. Detailed descriptions can be found in the notes to this R file.

2.	2-Sample Filtering-LUAD/LUSC/SKCM.R: This file is used to filter samples for the TPM data.

3.	3-Clinical Record-LUAD/LUSC/SKCM.R: This file is used to process clinical and exposure history variables for the selected samples.

4.	4-Gene Filtering and Pathways.R: This file links Ensembl ID and Entrez ID of genes and their involved pathways.

5.	5-Gene and Clinical merge.R: This file merges the TPM data and clinical records data and also generates the data frame format and matrix format data.

(Survival Analysis: survival analysis on processed TCGA data and selecting genes significantly impact survival)

6.	6-Marginal Survival Analysis.R: This file is used to perform marginal Cox proportional hazards model analysis.

7.	7-Joint Survival Analysis.R: This file is used to perform joint Cox proportional hazards model analysis with elastic net penalty.

8.	8-Robust Survival Analysis.R: This file is used to perform robust estimation using both least absolute deviations (LAD) and LAD with lasso penalty assuming an accelerated failure time (AFT) model.

(Replicability Analysis: replicability analysis for three measures and visualizing the generated measures)

9.	9-Replicability Analysis.R: This file includes codes for evaluating replicability of two selected gene sets based on three proposed replicability measures, as well as calculating measures under the null.

10.	10-Visualize Replicability.R: This file is used to generate plots and summary statistics of the replicability measures.

(Simulation: simulation studies to explore the properties of replicability measures, corresponds to manuscript section 5)

11.	11-Simulation (mimic true).R: This file includes simulation studies in section 5.1, which generates simulated data mimicking TCGA data and performs the replicability analysis on the results from the marginal and joint survival analysis.

12.	12-Simulation (whole vs true).R: This file includes simulation studies in section 5.2, which evaluates the impact of signal levels of genes on the replicability analysis results. This file also includes codes for plots.


