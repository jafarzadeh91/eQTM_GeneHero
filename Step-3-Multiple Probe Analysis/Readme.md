Readme
================
Tiam Heydari
2018-04-02

step3.Rmd and step3.md : collective effect of methylation probes on gene expression value
--------------
- These files contain the codes to the procedures in step 3.
- Use linear multiple regression model,
- Feature selection: Most significant probe, Forward probe selection, Backward probe selection (MASS), Lasso- (GLMnet), and Full probes.
-Use 80% of subjects as train and 20% as test. We train the model using 3-fold cross validation.

step3_files (Images Folder)
--------------

- The directory contains all images that are outputs from step3.md

step3-nonlinear.Rmd and step3-nonlinear.md : Non-linear models to predict gene expression level
--------------
- These files contain the codes related to the model extension to non-linear regression.
- Use neural network and conditional inference tree,

step3-nonlinear (Images Folder)
--------------
- The directory contains all images that are outputs from step3.md


Data
--------------

### rosmap_postprocV2.RData

- The data sets produced here consists of the results from Step 1. The .RData folder consists of three data sets: probes_subjects, subjects_genes, and probes_genes_distance. The probes_subjects data consists of DNA methylation probe data, subjects_genes consists of gene expression values, and probes_genes_distance consists of the distances on probes and genes where 0 represents a distance over 1Mb. Since it is too large to be uploaded in github, the data set could be found in this Google Drive Folder:

https://drive.google.com/drive/folders/1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY?usp=sharing

### cor_test_results_PCA_lapply_V4.rds

- The data set consists of the results on all correlation tests between probe and gene pair. This data set is used in Step 3 and Step 4. It is too large to be uploaded in github, so the link to the data set could be found here:

https://drive.google.com/drive/folders/1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY 

### probes_subjects_PCA_adjusted_V4.RDS

- The data set consists of probes_subjects from rosmap_postprocV2.RData which is the DNA Methylation probe data adjusted for 5 PCs. The data set could be found in the link here:

https://drive.google.com/drive/folders/1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY 

### subjects_genes_PCA_adjusted_V4.RDS

- The data set consists of subjects_genes from rosmap_postprocV2.RData which is the gene expression data adjusted for 30 PCs. The data set could be found in the link here:

https://drive.google.com/drive/folders/1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY 

