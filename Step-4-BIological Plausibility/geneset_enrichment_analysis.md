Step4- Biological Plausibility-Gene Set Enrichment Analysis
================
Zohreh Sharafian
March 27, 2018

Next, we study the relationship between numbers of significant methylated CpG sites corresponding to the genes and the function of those genes. In other words, we check if those genes that correspond to the multiple probes are multifunctional ane are more involved in biological process.

We load related library packages to understand and verify the biological aspects of our analysis.We use "ermineR" package for gene set enrichment analysis(Gillis et al., 2010).

### Loading data

WE load dataset from eQTM analysis that we have peroformed in Step 2 including all the adjusted P values for analyzing single probes.

``` r
 original_data <- readRDS("C:/Users/zohre/Desktop/Repo_team_Gene_Heroes/Data/cor_test_results_PCA_lapply_V4.rds")

head(original_data)
```

    ##                          gene      probe    estimate    pvalue
    ## cor  RAB4B:ENSG00000167578.11 cg25697727  0.04075525 0.3724587
    ## cor1 RAB4B:ENSG00000167578.11 cg02686662 -0.01476077 0.7467678
    ## cor2 RAB4B:ENSG00000167578.11 cg14319773  0.04823112 0.2911267
    ## cor3 RAB4B:ENSG00000167578.11 cg05498041 -0.05633210 0.2174921
    ## cor4 RAB4B:ENSG00000167578.11 cg14583103  0.06417159 0.1599679
    ## cor5 RAB4B:ENSG00000167578.11 cg18074151  0.03974841 0.3843954
    ##      adjusted.pvalue
    ## cor        0.9794159
    ## cor1       0.9957041
    ## cor2       0.9728284
    ## cor3       0.9633156
    ## cor4       0.9527597
    ## cor5       0.9801636

### Modifiying data

We remove the second part of gene names as "ermineR" accepts the HGNC symbol of the genes.

``` r
modified_data <- original_data

modified_data$new_gene_name <- sub(':.*$',"",modified_data$gene)

modified_data <- modified_data[c("gene","new_gene_name", "probe","estimate", "pvalue", "adjusted.pvalue")]

head(modified_data)
```

    ##                          gene new_gene_name      probe    estimate
    ## cor  RAB4B:ENSG00000167578.11         RAB4B cg25697727  0.04075525
    ## cor1 RAB4B:ENSG00000167578.11         RAB4B cg02686662 -0.01476077
    ## cor2 RAB4B:ENSG00000167578.11         RAB4B cg14319773  0.04823112
    ## cor3 RAB4B:ENSG00000167578.11         RAB4B cg05498041 -0.05633210
    ## cor4 RAB4B:ENSG00000167578.11         RAB4B cg14583103  0.06417159
    ## cor5 RAB4B:ENSG00000167578.11         RAB4B cg18074151  0.03974841
    ##         pvalue adjusted.pvalue
    ## cor  0.3724587       0.9794159
    ## cor1 0.7467678       0.9957041
    ## cor2 0.2911267       0.9728284
    ## cor3 0.2174921       0.9633156
    ## cor4 0.1599679       0.9527597
    ## cor5 0.3843954       0.9801636

### Preparing the inputs

We extract three columns including "new gane name", "probe", and "adjusted P-value" for all the genes in our dataset. Next, we count the number of probes corresponding to a particular gene and assign zero to the probees that are not significant and have P-value&lt;0.05.

``` r
extracted_data <- modified_data[,c(2, 3,6)]

filtered_adjusted_pvalue <- extracted_data[extracted_data$adjusted.pvalue<0.05,]
#filtered_adjusted_pvalue= extracted_data
all_genes_names=unique(modified_data$new_gene_name);
remaining_names=setdiff(all_genes_names,filtered_adjusted_pvalue$new_gene_name)
remaining_names=as.data.frame(remaining_names);
remaining_names$count=0;

filtered_adjusted_pvalue [1:20,]
```

    ##          new_gene_name      probe adjusted.pvalue
    ## cor63            RAB4B cg08800497    3.728588e-03
    ## cor260           RAB4B cg05033529    6.769418e-03
    ## cor297           RAB4B cg18434718    1.054311e-05
    ## cor523           RAB4B cg00524108    9.307602e-03
    ## cor740           RAB4B cg15744492    2.160051e-03
    ## cor790           RAB4B cg20176532    4.226004e-06
    ## cor8720         ZNF781 cg02715685    1.295830e-02
    ## cor18316        ZNF781 cg24088775    3.913203e-03
    ## cor2777         ZNF781 cg25612391    3.717254e-02
    ## cor3357         ZNF781 cg22628286    1.578355e-03
    ## cor4728         ZNF781 cg26539615    2.418363e-05
    ## cor5087         ZNF781 cg27068996    1.034188e-02
    ## cor5507         ZNF781 cg04607867    4.328839e-03
    ## cor5528         ZNF781 cg02844589    1.269695e-04
    ## cor2438           TH1L cg18369680    3.023043e-02
    ## cor9526   RP4-697K14.7 cg17069396    2.499666e-02
    ## cor23910  RP4-697K14.7 cg11024180    2.203888e-02
    ## cor9165   RP4-697K14.7 cg03307911    1.378800e-03
    ## cor10861  RP4-697K14.7 cg00152041    7.394562e-05
    ## cor52215        IL27RA cg23218533    3.970250e-04

We realize that there are 536 genes that have one or more significant probes.

``` r
nrow(filtered_adjusted_pvalue)
```

    ## [1] 1763

``` r
length(unique(filtered_adjusted_pvalue$new_gene_name))
```

    ## [1] 536

We generate a new dataset including the gene names and the counts of corresponding probes corresponding to each gene.Those genes that do not have any related probes get zero value.

``` r
library(tidyverse)

gene_list <- filtered_adjusted_pvalue %>% count(new_gene_name)
gene_list[1:15,]
```

    ## # A tibble: 15 x 2
    ##    new_gene_name     n
    ##    <chr>         <int>
    ##  1 A1BG              4
    ##  2 A4GALT            2
    ##  3 ABCA7             1
    ##  4 ABHD12            9
    ##  5 ABHD8             2
    ##  6 AC004258.1        1
    ##  7 AC004696.2        2
    ##  8 AC005003.1        1
    ##  9 AC006538.1        3
    ## 10 AC010335.1        2
    ## 11 AC012309.5        9
    ## 12 AC012615.1        9
    ## 13 AC016629.2       15
    ## 14 AC023490.1        1
    ## 15 AC092295.7        1

``` r
colnames(remaining_names)=c("new_gene_name","n")
gene_list=rbind(gene_list,remaining_names)
#gene_list$n=exp(gene_list$n)
class(gene_list)
```

    ## [1] "tbl_df"     "tbl"        "data.frame"

We provide what ermineJ requires as inputs including:

1.The gene ontology (GO) XML file 2.The annotation file The annotation file and GO XML files together provide the gene sets. 3.The list of genes and corresponding scores to be analyzed.

``` r
if (!file.exists("GO.xml")) { goToday("GO.xml") }

ermineInputGeneScores <- gene_list %>% 
  mutate(absolute_n = abs(n)) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_n)) %>% 
  column_to_rownames("new_gene_name")
```

    ## Warning: package 'bindrcpp' was built under R version 3.4.3

We use the Precision-Recall method to perform the enrichment analysis.

``` r
enrichmentResult <- precRecall(scores = ermineInputGeneScores, 
                               scoreColumn = 1, 
                               bigIsBetter = TRUE, 
                               annotation = "Generic_human", 
                               aspects = "B", 
                               iterations = 10000, 
                               geneSetDescription = "GO.xml") 
```

The result show the "correctedMFPvalue"=1 for all the genes regardless of the number of corresponding probes.This result shows that there is not any relationship between the number of significant probes corresponding to the genes and the multifunctionality of the genes as there is no significant correctedMFPvalue for the gene sets.

``` r
enrichmentResult$results %>% arrange(CorrectedPvalue) 
```

    ## # A tibble: 792 x 12
    ##    Name         ID     NumProbes NumGenes RawScore    Pval CorrectedPvalue
    ##    <chr>        <chr>      <int>    <int>    <dbl>   <dbl>           <dbl>
    ##  1 oxoacid met~ GO:00~        68       68   0.0784 2.20e-3           0.286
    ##  2 fatty acid ~ GO:00~        29       29   0.0622 2.90e-3           0.324
    ##  3 carboxylic ~ GO:00~        64       64   0.0783 2.20e-3           0.344
    ##  4 monocarboxy~ GO:00~        37       37   0.0695 2.00e-3           0.390
    ##  5 organic aci~ GO:00~        69       69   0.0783 4.60e-3           0.449
    ##  6 oxidation-r~ GO:00~        82       82   0.0921 6.00e-4           0.469
    ##  7 lipid metab~ GO:00~       111      111   0.105  2.00e-3           0.521
    ##  8 cellular li~ GO:00~        88       88   0.0922 1.40e-3           0.547
    ##  9 drug metabo~ GO:00~        52       52   0.0639 6.90e-3           0.599
    ## 10 response to~ GO:00~        23       23   0.0502 9.90e-3           0.703
    ## # ... with 782 more rows, and 5 more variables: MFPvalue <dbl>,
    ## #   CorrectedMFPvalue <dbl>, Multifunctionality <dbl>, `Same as` <chr>,
    ## #   GeneMembers <chr>
