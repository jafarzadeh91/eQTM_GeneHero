
**Studying the effect of DNA methylation on gene expression in human dorsolateral prefrontal cortex**
=====================================================================================================

This is the private repository of team **Gene Heroes** for STAT 540/ BIOF 540/ GSAT 540 term project. 


Members |  Department
---------|------------
Tiam Heydari | *Department of Biophysics*
Sina Jafarzadeh|	*Department of Bioinformatics*
Lisa Leung | *Department of Statisitcs*
Zohreh Sharafian|	*Department of Experimental Medicine*
Hiwot Tafessu |*Department of Statistics*


**Task Division Brief description of directories:**
---------------------------------------
**Step 0: Initial Idea** - Includes bringing up the initial idea, explaining different aspects of the question and the associated data, and introducing relevant papers for future research on the question by other members __Sina contributed to this step.__ 

**[Step 1: Data Processing](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-1-Data%20Processing)** - Includes description of data, codes for data preprocessing in preparation for data analysis and image outputs from this step. __Sina contributed in both the research and coding in Step 1.__ 

**[Step 2: Single Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-2-Single%20Probe%20Analysis)** - This directory includes code and images on PCA analysis, eQTM analysis and correlation results of each probe and gene pair. __Hiwot contributed in researching the topic and Lisa contributed in the coding.__

**[Step 3: Multiple Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-3-Multiple%20Probe%20Analysis)** - This directory will include the code and images on the exploration of multiple probe regression (full model using all probes that are significant), nonlinear regression, and LASSO variable selection, Stepwise (Forward and Backward) variable selection, neural network, and conditional inference tree method. __Tiam contributed researching and coding using regression modelling and variable selection methods for Step 3 (full model, nonlinear regression, Lasso and Stepwise methods), and Sina contributed in the researching and coding in neural network and  conditional inference tree method.__

**[Step 4: Biological Plausability](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-4-BIological%20Plausibility)** - This directory includes the visualizations and explanation behind biological reasonings and gene enrichment analysis. __For gene enrichment analysis and the analysis of spatial pattern of methylation probes, Sina contributed in the research and idea in addition to collaboration with Zohreh on coding part, and Zohreh contributed in the coding and writting aspect of the them. For analysis of the distribution of methylation probes in different chromatine states (ChromHMM), sina contributed to the research, coding and writing part.__

**[shinyApp](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/shinyApp)** - Data visiualization app to supplement analysis outputs and visualize expression and methylation data for selected genes. __Hiwot coded the shinyApp program.__ 

**[zz_Deliverable](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/zz_Deliverable)**:

- [Project proposal](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Group_Project_Proposal_Gene_Heroes.md) - __Lisa prepared the proposal, Hiwot and Tiam contributed to the ideas, then Zohreh and Sina finalized the proposal.__
- [Progress report](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Progressreport.md) - __Lisa prepared the report, Hiwot and Tiam contributed to the ideas, then Zohreh and Sina finalized the report.__
- [Poster presentation](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/zz_Deliverable/Final_Poster.pdf) - __Zohreh and Sina prepared and finalized the poster, then Hiwot and Tiam provided organized all the graphs plots.__

**Repo Preparation** - Organizing and preparing the entire repository at presentable level. __Hiwot and Lisa organized the Repo. Zohreh and Sina wrote Introduction, Background and Data Description. Lisa prepared the summary analyses but are finalized by everyone accordingly to the steps they have contributed. Zohreh wrote Bibliography.__

**Introduction**
----------------
DNA methylation is an epigenetic mechanism which plays a role in regulating tissue specific gene expression(Acharya, Owzar, & Allen, 2017). DNA methylation mostly occurs when a cytosine base is next to a guanine base, forming a CpG site. Early studies showed that methylation of CpG sites prevent the expression of genes. However, recent studies revealed that methylation can be linked with both decreasing and increasing of gene expressions(Gutierrez-Arcelus et al., 2013). 
In this study, we characterize the relationship between DNA methylation and gene expression in a quantitative manner, i.e. evaluating the statistical significance of the effect of methylated CpG sites on gene expressions colloquiality known as **expression Quantitative Trait Methylation (eQTM)** in human dorsolateral prefrontal cortex (DLPFC).


**Goal**
----------------

There are a few studies quantizing the relationship of methylation and gene expression in a tissue specific manner(Bonder et al., 2017). In this project, we perform eQTM  of two data sets derived from the dorsolateral prefrontal cortex (DLPFC)(Bennett, Schneider, Buchman, et al., 2012; Bennett, Schneider, Arvanitakis, & Wilson, 2012). 
**1)** We analyze the correlation of each gene and the probes located near that gene. **2)** We extend our model to include the collective effect of multiple CpG sites on the expression of genes. **3)** We quantitatively assess the roles of chromatin states and CpG site distance from gene on the regulatory relationship between methylation probes and a gene expression on brain tissue. This provides us more biological intuition behind the regulatory processes in human brain.



![](https://raw.githubusercontent.com/STAT540-UBC/Repo_team_Gene_Heroes/master/zz_Deliverable/figone.PNG?token=AVAzvOigxIEpdo_qYdVL1Q2YBvuf6rGjks5aznpIwA%3D%3D)
**Figure 1. Graphical overview of our data analysis pipeline.** We first correct the data for confounders and outliers. Next, we apply eQTM analysis to assess significant correlations between the pairs of methylation probe and gene. We then extend the analysis to multivariate linear and non-linear regression models. We verify the biological plausibility of results using other data types including probe-gene distances and chromatin states.


**The Dataset**
----------------
We use a combination of two publicly available data sets, ROSMAP, namely as The Religious Orders Study (ROS)(Bennett, Schneider, Arvanitakis, et al., 2012) and The Memory and Aging Project (MAP)(Bennett, Schneider, Buchman, et al., 2012). Subjects in ROSMAP were collected in 1994 and 1997 respectively by the same researchers. Both data sets are longitudinal cohort studies containing DNA methylation probes and gene expression data, and phenotypes of diseases related to aging. However, we do not use the phenotypes for the purposes of this project due to the confidentiality agreement. The gene expression data is RNA-sequencing data extracted from dorsolateral prefrontal cortex from 540 subjects using Illumina HiSeq with 101-bp paired-end reads. The DNA methylation data was generated using 450K Illumina array from 740 subjects. There are total of 468 samples individuals who have both the data. The data is Quality Controlled and corrected for common batch effect and cofounders(Bennett, Schneider, Buchman, et al., 2012; Bennett, Schneider, Arvanitakis, et al., 2012). Due to computational limitation, we study chromosome 18-22 out of the whole dataset.


**Summary of analysis and major results**
------------------------

**[Step 1: Data Processing](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1.md)**: ROSMAP raw data is preprocessed to prepare for the single probe and multiple probe analyses. We started with a gene expression data set, DNA methylation data set and its corresponding CpG sites. For this particular assignment, we are only interested in the probes within a 1Mb distance around the gene of interest. 

In order to demonstrate our multiple probe analysis, we have to handle with the issue of multicollinearity between probes.  Hence, we decide to pick one probe out of a cluster of probes that are strongly correlated with each other (ie. probes that have strong linkage disequilibrium). We used A-clust algorithm to solve this issue. Based on the resulting plots of [number of non-singletone clusters](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-1.png) and [average intra-cluster correlation](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-1-Data%20Processing/Step1_files/figure-markdown_github/unnamed-chunk-1-2.png), we decide that correlation threshold of 0.8 and a base-pair distance of 1000 are the ideal parameters to separate the probes. In return, we are left with three data sets for Step 2 and 3: [probes_genes_distance, probes_subjects and subjects_genes](https://drive.google.com/file/d/1Ieze9KwSy5UL9c6Vt5uFF9ta4U9kvkoe/).

**[Step 2: Single Probe Analysis](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2.md)**: Using the [data processed in Step 1](https://drive.google.com/file/d/1Ieze9KwSy5UL9c6Vt5uFF9ta4U9kvkoe/) , the PCA is used to adjust for any hidden confounders that could exist in our data sets (in both gene expression and DNA methylation). To choose the number of PCs adjusted for both gene expression data and DNA methylation probe data, we plotted the [heatmap](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/prop_of_genes_sig_on_50_genes_V4_fdr.png) on the proportion of genes that have at least one significant probe with FDR < 0.1 based on number of PCs adjusted for both data sets(Ma & Dai, 2011). Note that the heatmap is based on 50 genes. For exploration, we also plotted the percentage of cumulative variance that additional PCs are adjusted for in both [gene expression data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/cum_var_explained_50_V4_gdata.png) and [DNA methylation data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/cum_var_explained_50_V4_probe.png). The optimal number of PCs adjusted will be the highest number of significant values based on the heatmap. 

As a [result](https://drive.google.com/open?id=1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY), we concluded that the optimal DNA methylation probe data set should be adjusted for five PCs and gene expression data set is adjusted for 30 PCs. To observe [the effects of adjusting for number of PCs in the gene expression data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/gene_C21orf56_gdata.png), we took a look at a gene that show the most number of significant probes. Similarly, [the effects of adjusting for number of PCs in DNA methylation data](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-2-Single%20Probe%20Analysis/Step-2_files/gene_C21orf56_probes.png) are explored in the same gene.  


**[Step 3: Multiple Probe Analysis on Statistical Models](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis)**: Based on the single probe analysis, we would like to explore what happens when we modeled linear regression with more than one probe for each gene expression. To achieve this, we first divided our data set to 80% train and 20% test. Aslo we used 3 fold cross validation during the analysis. Then we ran a single linear regression (using the probe with the smallest FDR value) and compared with the multiple linear regression using all significant probes (FDR < 0.05). In addition, we also ran Forward and Backward stepwise variable selection method for each gene, LASSO variable selection method. to assess the effect of higher degree fitting on our data, we also checked the higher degree polynomial fits on the single probes.

The [results shows that the full probes set](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3_files/figure-markdown_github/unnamed-chunk-15-1.png) is the best model for more than 50% of the genes. interestingly the second best used the most significant probes transformed by a high degree function.

To deal with the danger of overfitting we also [plot the ratio of test to train error](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3_files/figure-markdown_github/unnamed-chunk-16-6.png) and we found that it is symmetrically distributed. around 1 which indicate that we didn't overfit the data.

The results of this part suggest that firstly, utilizing the full set of probes provides the least value of errors in comparison the picking a subset of probes as covariates. Secondly, the non-linear transformation of the most significant CpG probe provides less error in comparison the linear model by that CpG probe. For this purpose, [we chose two non-linear model to consider the non-linearity of the relationship between inputs and outputs](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3-nonlinear.md). We train a neural network with one hidden layer and at most 5 neurons (due to the computational limitations) and we use the logistic function as activation function to involve non-linearity. [Here](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3-nonlinear_files/figure-markdown_github/unnamed-chunk-5-1.png), we see the error of neural network model on different settings (number of neurons and regularization function). We see that the model is robust against overfitting and increasing the number of neurons improves the error. A similar diagram for conditional inference tree [here](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3-nonlinear_files/figure-markdown_github/unnamed-chunk-5-2.png) shows the same trend about the model robustness against overfitting. Consistent with observation in linear models, we see that including all significant CpG probes for a gene (instead of just most significant CpG probe) leads to improvement in gene expression estimation. [here](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-3-Multiple%20Probe%20Analysis/step3-nonlinear_files/figure-markdown_github/unnamed-chunk-6-1.png), we plot the contribution of all two significant probes for a random gene in different models showing the contribution of each probe to the overall r-square. We finally compare the mse error of models on test data. While The average error of neural network and conditional inference tree are close (i.e. 0.033 and 0.35) they do not significantly outperform the linear models as we expected initially. In fact, the average mse for genes across their best-fitted model is 0.0358 which is just slightly higher than our neural network. We expect the improvement increases if we increase the number of hidden layers and neurons in the neural network and we are going to address it in future works.

**[Step 4: Biological Plausability](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/tree/master/Step-4-BIological%20Plausibility)**: Based on the single probe analysis in Step 1, we assess [the spatial distribution of significant methylation probes](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4-spatial_distribution_of_methylated_CpG_sites.md) according to their distance from genes at FDR < 0.05. The [histogram](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4-spatial_distribution_of_methylated_CpG_sites_files/figure-markdown_github/unnamed-chunk-6-1.png) shows that the majority of methylated CpG sites are concentrated near the transcription start sites (TSS) of genes. The [bar graph](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4-spatial_distribution_of_methylated_CpG_sites_files/figure-markdown_github/unnamed-chunk-9-1.png) shows that the methylated CpG sites with significant P-values are located close to TSS. This result is consistent with previous data which showed that CpG sites are relatively enriched around TSS of genes (Numata et al., 2012; Numata, Ye, Herman, & Lipska, 2014; Saxonov et al., 2006).
We also analyze the [distribution of significant methylation probes in different chromatin states](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4-chromatine_state_pattern.md) using [ChromHMM](http://compbio.mit.edu/ChromHMM/) (Ernst & Kellis, 2012). [The bar graph](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/step4-chromatine_state_pattern_files/figure-markdown_github/unnamed-chunk-1-1.png) show the frequency of total number of positive and negative regulating CpG sites for different chromatin states. It shows that TSS enriches with the highest number of significant methylation probes compared to other chromatin states.
Moreover, we see that the number of negative regulating methylation probes in TSS region is higher compared to the positive ones. Biological experiments also verify that most of the methylated CpG sites near the TSS and gene promotor
repress the expression value of the gene resulting in a down-regulation effect(Numata et al., 2012; Numata, Ye, Herman, & Lipska, 2014; Saxonov et al., 2006).
We also perform [“Gene Set enrichment analysis”](https://github.com/STAT540-UBC/Repo_team_Gene_Heroes/blob/master/Step-4-BIological%20Plausibility/geneset_enrichment_analysis.md) using [erminR](https://github.com/STAT540-UBC/STAT540-UBC.github.io/blob/master/seminars/seminars_winter_2017/seminar10/sm10_gene_set_analysis.md) to find the relationship between the number of significant probes and the complexity of gene regulation using erminR package (Gillis et al., 2010). The analysis shows that there is not
any significant relationship between the numbers of significant methylated CpG sites corresponding to genes and the multifunctionality.


**Discussion**
----------------
Our study presents a large-scale analysis on the regulatory effect of methylation CpG sites on the expression value of genes in brain dorsolateral prefrontal cortex. We develop a set of regression models to assess the relationship between expression value of each gene and the methylation value of CpG sites. We investigate the distribution pattern of CpG sites considering the chromatin state and CpG-Gene distance. these findings are consistent with previous biological studies(Gutierrez-Arcelus et al., 2013; Numata et al., 2012; Saxonov, Berg, & Brutlag, 2006). 
We observed that despite our initial expectation, the non-linear models (NN and CIT) did not provide better estimations in comparison to linear regression. We suspect that using more complicated models, e.g. increasing the number of hidden layers in NN may capture the potential non-linearity between the CpG-gene relationship. A line of research that we will address in future.



**References**
----------------

Acharya, C. R., Owzar, K., & Allen, A. S. (2017). Mapping eQTL by leveraging multiple tissues and DNA methylation. BMC Bioinformatics, 18(1), 1–20.

Bennett, D. A., Schneider, J. A., Arvanitakis, Z., & Wilson, R. S. (2012). overviw and findings from the religious orders study. Current Alzheimer Research, 29(10), 628–645.

Bennett, D. A., Schneider, J. A., Buchman, A. S., Barnes, L. L., Boyle, P. A., & Wilson, R. S. (2012). Overview and findings from the rush Memory and Aging Project. Current Alzheimer Research, 9(6), 646–63.

Bonder, M. J., Luijk, R., Zhernakova, D. V., Moed, M., Deelen, P., Vermaat, M.,Heijmans, B. T. (2017). Disease variants alter transcription factor levels and methylation of their binding sites. Nature Genetics, 49(1), 131–138.

Ernst, J., & Kellis, M. (2012). ChromHMM: Automating chromatin-state discovery and characterization. Nature Methods, 9(3), 215–216.

Gillis, J., Mistry, M., & Pavlidis, P. (2010). Gene function analysis in complex data sets using ermineJ. Nature Protocols, 5(6), 1148–1159.

Gutierrez-Arcelus, M., Lappalainen, T., Montgomery, S. B., Buil, A., Ongen, H., Yurovsky, A., Dermitzakis, E. T. (2013). Passive and active DNA methylation and the interplay with genetic variation in gene regulation. eLife, 2013(2), 1–18.

Ma, S., & Dai, Y. (2011). Principal component analysis based Methods in bioinformatics studies. Briefings in Bioinformatics, 12(6), 714–722.

Numata, S., Ye, T., Herman, M., & Lipska, B. K. (2014). DNA methylation changes in the postmortem dorsolateral prefrontal cortex of patients with schizophrenia. Frontiers in Genetics, 5(JUL), 1–7.

Numata, S., Ye, T., Hyde, T. M., Guitart-Navarro, X., Tao, R., Wininger, M.,Lipska, B. K. (2012). DNA methylation signatures in development and aging of the human prefrontal cortex. American Journal of Human Genetics, 90(2), 260–272. 

Saxonov, S., Berg, P., & Brutlag, D. L. (2006). A genome-wide analysis of CpG dinucleotides in the human genome distinguishes two distinct classes of promoters. Proceedings of the National Academy of Sciences, 103(5), 1412–1417.


---------------------------------------

........

