Progress Report
================
Zohreh Sharafian
March 19,2018

### What has changed based on the final proposal (1 pt.)

-   **Did your dataset change? If so, why?**

The dataset we planned to use have not changed. However, we used a subset (chromosomes 18 - 22) of the data set due to the computational resource limitation. In step 3 and 4, we might also decide to narrow down analysis to the specific genes to meet the available resources.

-   **Have you decided to do a different analysis than what was mentioned in your proposal? If so, Why?**

We decided to keep the same analysis mentioned in our proposal.

-   **Are there any changes in task assignments of group members?**

There have been slight changes in task assignments due to feasibility. As of present, Sina has finished Step 1 of the proposal; Lisa and Hiwot have finished Step 2 of the proposal. We are currently on Step 3 of the proposal which is assigned to Zohreh and Tiam.

### What is the progress of the analyses (4 pts.)

**1.1 Preprocessing**

Our data consists of the gene expression and methylation measurements for 508 and 702 subjects, respectively. Apparently, it imposes an unnecessary huge number of hypothesis to test the correlation between each single CpG site and gene pair. So, We solely analyze the probes located in a 1MB window around the gene of interest.

**1.2 Multicollinearity of Methylation Probes**

It is known that methylation probes can be in strong Linkage Disequilibrium (LD) with their neighbouring probes. So, we need to reduce the number of probes to a set of unrelated ones. It helps us to decrease the computational cost of the model. Moreover, It is necessary to have uncorrelated probes in the third phase of the project where we provide a linear predictor of multiple probes for the gene of interest. There are a few tools developed for this purpose (1,2). In this project, we use A-clust algorithm provided by Harvard University (1).

This algorithm first merges each pair of nearby probes and all of the probes located between them if the correlation between that pair is more than corr\_dist\_thresh. Nearby pairs are pairs of probes with a base pair distance less than *bp\_dist\_thresh*. In the next step, it merges each two neighbour probes if their correlation distance is less than *corr\_dist\_thresh*. We find the best values for these parameters by creating train and validation datasets. We choose the parameters based on train data and then check if they work on another validation data. We choose CpG sites on chromosomes 1-18 as train data and CpG sites in chromosome 19-22 as validation. We select the parameters that decrease the number of probes (i.e. create as much cluster as possible) and yield the highest intra-cluster correlation between probes.

![dsfsff](https://raw.githubusercontent.com/STAT540-UBC/Repo_team_Gene_Heroes/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-1.png?token=AVAzvNXCSrOSnYHW_KIwVEU1aw5ZMaK2ks5azrA4wA%3D%3D)

![dsfsff](https://raw.githubusercontent.com/STAT540-UBC/Repo_team_Gene_Heroes/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-2.png?token=AVAzvNb7Gn5h8iZjZa-RcOPAd0Sx9Wtwks5azrBGwA%3D%3D)

The plots suggest that *corr\_threshold* = 0.8 and *bp\_dist\_threshold* = 1000 are the best parameters in terms of the intra-cluster correlation/number of clusters trade-off. The authors of A-clust paper also suggested this parameter in their experiment on human methylation data, supporting the chosen parameters in our problem.

**1.3 Data Conversion**

In this part, we convert the raw files to R language data format to make it easier to handle. We also combine the raw files to organize the data. In order to lower the computational cost, we limit the analysis to chromosomes 19-22. Below is the list of final data structure saves together as rosmap\_postdocV2.RData file:

-   **probes\_genes\_distance**: a sparse matrix created by Matrix library showing the distance of each cpg probe to each gene in its neighbourhood.

-   **probes\_subjects**: a data frame showing the Z-scored value for each cpg site across all subjects.

-   **subjects\_genes**: a data frame showing the expression value of each gene across all subjects.

**2.1 Single Probe analysis**

After the preprocessing methods in Step 1, we kept the original plan in Step 2 for performing single probe analysis. However, in addition to the QC steps performed in Step 1, we needed to address the issue of hidden confounding factors before moving forward to Step 2. So, we further explored our data by plotting the estimated proportion of genes that have at least one significant eQTMs corrected by FDR p-value at 0.05 threshold and also plotted the proportion of variance that will be explained. Based on our observations and asking the class instructor, we decided that adjusting for 3 PCs that explained approximately 30% of the variation was sufficient and accordingly adjusted our data.

The result shows the proportion of genes with at least one significant eQTM by looking at FDR &lt; 0.05 against the number of PCs adjusted in the gene expression data set. Due to computational timing issues, we are only able to look at a set of 50 random genes to plot instead of looking at all the genes. This plot is slightly difficult to interpret, but it seems like the proportion of genes increase in a general direction when more PCs are adjusted.

![](https://raw.githubusercontent.com/STAT540-UBC/Repo_team_Gene_Heroes/ecf9d6c878e201e629133ec8a106cc8325f3ff92/images/cum_var_explained_all_V2.jpg?token=AVAzvOMI69ZW0qMEFdEqzau2XsAotlvFks5aztCKwA%3D%3D)

The figure above shows the cumulative variance explained against the number of principal components adjusted. We decide to adjust for three principal components with approximately 0.3 variance explained overall.

-   **What R packages or other tools are you using for your analyses? You do not need to provide your scripts in your report.Provide the links to any markdown reports within your repo to refer to the relevant analysis.**

Being the raw data in Matlab and Excel formats, we used various data conversion base libraries to transform data to appropriate R format. In order to represent sparse matrices, we used Matrix library. We adjusted A-clust external R package to address the multicolearity issue in our data as well. We used various statisticall methods from R base libraries to perform principle component analysis, hypothesis testing, etc. Aditional information is availabe in relevant Markdown files, [step1](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1.md) and [step2](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2.md).

-   **Provide references.**

1.  Sofer, T., Schifano, E. D., Hoppin, J. A., Hou, L. & Baccarelli, A. A. A-clustering: A novel method for the detection of co-regulated methylation regions, and regions associated with exposure. Bioinformatics 29, 2884-2891 (2013).
2.  Haque, M. M., Nilsson, E. E., Holder, L. B. & Skinner, M. K. Genomic Clustering of differential DNA methylated regions (epimutations) associated with the epigenetic transgenerational inheritance of disease and phenotypic variation. BMC Genomics 17, 1-13 (2016).

### Results (2 pts.)

-   **What are your primary results? Were you able to answer your hypothesis? Did you have any positive results? If no, postulate a discussion as to why that may be. Provide plots and/or tables to present your results.**

Our results from the single probe analysis can help us identify which methylation probes are significantly associated with expression levels. In order to evaluate our findings from this analysis, we perform literature review for selected genes and verifying the biological plausibility. As an example, we have included the table below to show the result for a randomly chosen gene 'RAB4B' and the top 10 significantly associated methylation probes.

![](https://raw.githubusercontent.com/STAT540-UBC/Repo_team_Gene_Heroes/ecf9d6c878e201e629133ec8a106cc8325f3ff92/images/tableresult.PNG?token=AVAzvMAiFdo3eOuUUEy7tIR41jVYcm4fks5aztBSwA%3D%3D)

In addition, we are currently working on multiple probe analysis to assess how combinations of methylation probes can better explain gene expression levels compared to single probe analysis.

-   **List some challenges that you encountered? How will you address them?**

The aims and the methodologies we planned on our proposal remain the same, and we have described in previous sections to address our objectives. However, a challenge we might face is biologically validating our findings from multiple probe analysis. In fact, we are looking for an appropriate tool or R library to study the biological relevance of the results.
