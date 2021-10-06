Group Project Proposal\_Gene Heroes
================
Zohreh Sharafian
February 14, 2018

Motivation and Background:
--------------------------

Addressing the relation of molecular level variations and complex diseases has been an interest of scientific community in the last decade. Genome and epigenome wide association studies have found tens of thousands of potential variations relevant to complex diseases in humans; However, most of them couldnot be exactly replicated due to the complications of human genome \[1,2\]. It is well known that DNA methylation affects the regulatory process of cell machinery, and subsequently the gene expressions \[3\]. However, very few studies investigate the correlation between DNA methylation and gene expression in a quantitative manner (eQTM), i.e. finding the CpG sites affecting the gene expression and the significance of their effects. Previous eQTM studies consider one CpG site per gene at a time. However, in reality, there might be a set of CpG sites that collectively explain the level of gene expression in a linear or nonlinear manner, while each site contributing a small amount to the overall effect.

**Objective:** In this project, we investigate which methylation probes are good predictors for the level of gene expression of different genes in human dorsolateral prefrontal cortex (DLPFC). First, we study one methylation probe at a time as an explanatory variable for each gene of interest. Then, we expand the model to leverage a set of methylation sites, instead of one CpG site, for each gene.

Division of Labour:
-------------------

-   Lisa Leung: R-Shiny/R code/Writing
-   Zohreh Sharafian: Background and Objective/R code/Poster/Writing
-   Sina Jafarzadeh: Background and Objective/R code/Poster/Writing
-   Tiam Heydari: Research/R code/Poster/Writing
-   Hiwot Tafessu: R-Shiny/R code/Writing

Dataset:
--------

We use a combination of two publicly available data sets, ROSMAP, namely as The Religious Orders Study (ROS) \[4\] and The Memory and Aging Project (MAP) \[5\]. Subjects in ROSMAP were collected in 1994 and 1997 respectively from the same researchers. Both data sets are longitudinal cohort studies containing DNA methylation probes and gene expression data, and phenotypes of diseases related to aging. However, we do not use the phenotypes for the purposes of this project due to the confidentiality agreement. The gene expression data is RNA-seq data extracted from dorsolateral prefrontal cortex from 540 subjects using Illumina HiSeq with 101-bp paired-end reads. The DNA methylation data was generated using 450K Illumina array from 740 subjects. There are a total of 468 samples individuals who have both the data1. The data is Quality controlled and corrected for common batch effect and co-founders \[4,5\].

Aims and Methodology:
---------------------

1.  **Preprocessing:** We will attempt to improve the quality of data by removing hidden confounders \[6\] using Principal Component Analysis (PCA) \[7\]. In addition, we will normalize the data and remove outliers \[8\]. Then, we reduce the number of correlated probes to avoid multicollinearity problem using the state-of-the-art open source tools such as A-clust \[9\] and the method by Haque *et al.* \[10\].
2.  **Single Probe analysis:** We use eQTM analysis to find significant effect between methylation probes and gene expression. We consider methylation probes within a 1 MB window. We correct the multiple hypothesis testing on results using the Family-Wise Error Rate (FWER) or False Discovery Rate techniques.
3.  **Multiple Probe analysis:** Based on our significant findings from the single probe analysis, we extend the model to a linear combination of methylation probes onto a gene (i.e. multiple linear regression). If time permits, we would explore other possible models such as nonlinear regression and neural networking methods aiming to consider the non-linear relationships
4.  **Biological Plausibility:** From the findings in the Multiple Probe analysis, we will pick a couple of genes and DNA methylation relationships and compare to previous findings in the literature. Through these reviews, we can understand and verify the biological aspect of our statistical findings.
5.  **Presentation:** Finally, we will develop an R-Shiny App and a poster as a report and a presentation.

References:
-----------

1.  Ng, B. et al. An xQTL map integrates the genetic architecture of the human brain-s transcriptome and epigenome. Nat. Neurosci. 20, 1418-1426 (2017).
2.  Feinberg, A. P. The epigenetic basis of common human disease. Trans. Am. Clin. Climatol. Assoc. 124, 84-93 (2013).
3.  Jaenisch, R. & Bird, A. Epigenetic regulation of gene expression: How the genome integrates intrinsic and environmental signals. Nat. Genet. 33, 245-254 (2003).
4.  Bennett, D. A., Schneider, J. A., Arvanitakis, Z. & Wilson, R. S. overviw and findings from the religious orders study. Curr. Alzheimer Res. 29, 628-645 (2012).
5.  Bennett, D. A. et al. Overview and findings from the rush Memory and Aging Project. Curr. Alzheimer Res. 9, 646-63 (2012).
6.  Listgarten, J., Kadie, C., Schadt, E. E. & Heckerman, D. Correction for hidden confounders in the genetic analysis of gene expression. Proc. Natl. Acad. Sci. 107, 16465-16470 (2010).
7.  Abraham, G. & Inouye, M. Fast principal component analysis of large-scale genome-wide data. PLoS One 9, 1-5 (2014).
8.  Barghash, A. & Arslan, T. Robust Detection of Outlier Samples and Genes in Expression Datasets. J. Proteomics Bioinform. 9, 38-48 (2016).
9.  Sofer, T., Schifano, E. D., Hoppin, J. A., Hou, L. & Baccarelli, A. A. A-clustering: A novel method for the detection of co-regulated methylation regions, and regions associated with exposure. Bioinformatics 29, 2884-2891 (2013).
10. Haque, M. M., Nilsson, E. E., Holder, L. B. & Skinner, M. K. Genomic Clustering of differential DNA methylated regions (epimutations) associated with the epigenetic transgenerational inheritance of disease and phenotypic variation. BMC Genomics 17, 1-13 (2016).
