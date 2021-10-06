data\_description
================
Sina Jafarzadeh
2018-02-14

Data Description
----------------

We will use a combination of two publicly available and preprocessed data sets, **ROSMAP**, namely as The **Religious Orders Study (ROS)** and **The Memory and Aging Project (MAP)**. Subjects in ROSMAP were gathered in 1994 and 1997 respectively from the same researchers. Both data sets are longitudinal cohort studies containing **DNA methylation probes** and **gene expression** data, and phenotypes of diseases related to aging. However, we do not use the phenotype information in this project due to the confidentiality aggreement. The gene expression data is RNA-seq data collected from dorsolateral prefrontal cortex tissue of 540 subjects using Illumina HiSeq technology with 101-bp paired-end reads. The DNA methylation data was generated using 450K Illumina array from 740 subjects. There are a total of 468 samples individuals who have both the data. The data is QC’ed and corrected for common batch effect and co-founders.

Data Structure
--------------

Below, we describe the data structures that are forming all the ROSMAP dataset.

**expressionAndPhenotype.mat** : This Matlab data file contains the gene expression data table for 504 subjects (rows) across 13484 genes (columns). This file also contains a 504 subjects by 13 phenotypes data table, showing the different phenotypes (e.g. smoking status, gender, etc) for all subjects.

**methylationSNMnorm.mat** : This Matlab data file contains a data table, showing the value measured for 420103 CpG sites (rows) for 702 subjects (columns).
