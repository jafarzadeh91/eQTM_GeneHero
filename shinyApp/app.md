shinyApp
================
Hiwot Tafessu

``` r
library(shiny)
library(shinythemes)
library(tidyverse) 
```

    ## -- Attaching packages ---------------------------------- tidyverse 1.2.1 --

    ## v ggplot2 2.2.1     v purrr   0.2.4
    ## v tibble  1.4.2     v dplyr   0.7.4
    ## v tidyr   0.8.0     v stringr 1.3.0
    ## v readr   1.1.1     v forcats 0.3.0

    ## -- Conflicts ------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(reshape2) 
```

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

``` r
library (corrplot)
```

    ## corrplot 0.84 loaded

``` r
library(psych)
```

    ## 
    ## Attaching package: 'psych'

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     %+%, alpha

``` r
#Local Directory (Please see Data file (Step-2))
subjects_genes<-as.data.frame(readRDS("~/Stat_540/zz_tafesu-hiwot_STAT540_2018/Group_work/subjects_genes_PCA_adjusted_V4.RDS"))
probes_subjects<-as.data.frame(readRDS("~/Stat_540/zz_tafesu-hiwot_STAT540_2018/Group_work/probes_subjects_PCA_adjusted_V4.RDS"))
cor_test_results<-readRDS("~/Stat_540/zz_tafesu-hiwot_STAT540_2018/Group_work/cor_test_results_PCA_lapply_V4.RDS")


#prep data for plot 
dat<- subjects_genes%>% 
  rownames_to_column("subject") %>% 
  melt(id = "subject")%>% 
  select(subject,
         gene = variable, 
         expression = value) 
dat$gene_name<- sub(':.*$',"",dat$gene)#shorten gene name from raw data

dat_cor<- cor_test_results
dat_cor$gene_name<-sub(':.*$',"",dat_cor$gene)#shorten gene name from raw data

gene.name  <- names(subjects_genes)
select10<- 1:10  
# function to id outlier subjects in plots
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(na.omit(x)) | x > quantile(x, 0.75) + 1.5 * IQR(na.omit(x)))
}


#-----------------------------------------------------------------
#From step-3
#plot the correlation between probes, genes 
plot_correlation <- function(designated.dataframe.for.specific.gene,number.of.picked.probes){
  plot(designated.dataframe.for.specific.gene[,1:min(number.of.picked.probes,5)], pch=1,main=gene.name)
  datacor = cor(designated.dataframe.for.specific.gene[1:number.of.picked.probes])
  corrplot(datacor, method = "color",addCoef.col="grey",number.cex= 7/ncol(designated.dataframe.for.specific.gene))
  auto.sel <- designated.dataframe.for.specific.gene[,1:number.of.picked.probes]
  pairs.panels(auto.sel, col="red")
}
#-----------------------------------------------------------------

#Prepare blank panels, input and output objects 
ui <- fluidPage(theme = shinytheme("slate"),
               
                tabPanel(titile=NULL,   titlePanel("ROSMAP Data (Chr18-22)"),
                         mainPanel(navbarPage(title = NULL,
                                      tabPanel("Homepage",  #brief description of shinyApp
                                                       h4("Team: Gene Heroes"),
                                                       h5("This Shiny app is prepared to supplement the term project for STAT 540 titled [Studying the regulatory effect of DNA methylation on gene expression in human dorsolateral prefrontal cortex (DLPFC)].  ROSMAP, is a dataset that resulted from a combination of two cohort studies, namely The Religious Orders Study (ROS) and The Memory and Aging Project (MAP), that contain gene expression and DNA methylation probe data extracted from the dorsolateral prefrontal cortex of aging subjects. This shiny app helps to visualize this data from 481 subjects, 1795 genes and 42813 probes. Please refer to https://github.com/STAT540-UBC/Repo_team_Gene_Heroes for more detail.")),
                                              
                                              tabPanel("Comparing Gene Expression",
                                                       plotOutput("per_gene_plot"),
                                                       selectInput( "select","Select other genes to compare? [Max=10]",c("Not Selected","Select other genes to compare")),
                                                       
                                                       conditionalPanel(condition = "input.select=='Select other genes to compare'", #first select random genes for display 
                                                                        selectInput("geneInput1", "1st Gene", c("None",levels(factor(dat$gene_name))),selected = 'AAR2'),
                                                                        selectInput("geneInput2", "2nd Gene", c("None",levels(factor(dat$gene_name))),selected = 'BACE2'),
                                                                        selectInput("geneInput3", "3rd Gene", c("None",levels(factor(dat$gene_name))),selected = 'CACTIN'),
                                                                        selectInput("geneInput4", "4th Gene", c("None",levels(factor(dat$gene_name))),selected = 'DLGAP4'),
                                                                        selectInput("geneInput5", "5th Gene", c("None",levels(factor(dat$gene_name))),selected = 'EID2'),
                                                                        selectInput("geneInput6", "6th Gene", c("None",levels(factor(dat$gene_name))),selected = 'ELAVL1'),
                                                                        selectInput("geneInput7", "7th Gene", c("None",levels(factor(dat$gene_name))),selected = 'FXYD5'),
                                                                        selectInput("geneInput8", "8th Gene", c("None",levels(factor(dat$gene_name))),selected = 'HRH3'),
                                                                        selectInput("geneInput9", "9th Gene", c("None",levels(factor(dat$gene_name))),selected = 'LPAR2'),
                                                                        selectInput("geneInput10", "10th Gene", c("None",levels(factor(dat$gene_name))),selected = 'MYO9B')
                                                       )
                                              ),
                                              tabPanel("eQTM", 
                                                       h4("The Spearman rank-order correlation"),
                                                       
                                                       selectInput("geneInput_cor", "Please select gene:", c("None",levels(factor(dat$gene_name)))),
                                                       tableOutput("corTable")
                                                       
                                                       
                                              ),
                                              tabPanel("Multiple Probes", 
                                                       
                                                       h5("Correlation Plots:"),
                                                       
                                                       selectInput("geneInput_cor2", "Please select gene:", c("None",levels(factor(gene.name )))),
                                                       
                                                       sliderInput("topProbes", "Number of probes:",
                                                                   min = 0, max = 50, value = 5),
                                                       plotOutput("corPlot"),
                                                       plotOutput("corPlot2")
                                                       
                                                       
                                                       
                                              )
                         )         
                         )
                         
                )
)



server <- function(input, output) {
  #output on Comparing Gene Expression panel
  output$per_gene_plot <- renderPlot({
    per_gene_data <-
      dat %>%
      filter(  gene_name == input$geneInput1[1] | gene_name == input$geneInput2[1] | gene_name == input$geneInput3[1]| 
                 gene_name == input$geneInput4[1] | gene_name == input$geneInput5[1] | gene_name == input$geneInput6[1]|
                 gene_name == input$geneInput7[1] | gene_name == input$geneInput8[1] | gene_name == input$geneInput9[1]|
                 gene_name == input$geneInput10[1]) %>% 
      mutate(is.outlier=ifelse(is.na(expression), FALSE, is_outlier(na.omit(expression)))
      ) %>% 
      mutate(outlier=ifelse(is.outlier == T, subject, as.numeric(NA)))
    
    ggplot( per_gene_data, aes(y=expression, x=gene_name)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 1)+
      ggtitle("Normalized distribution of gene expression for selected genes") +
      labs(x = "Selected Gene", y="Expression")+
      theme(plot.background = element_rect(fill = 'lightcyan4', colour = 'black'))
    
  })
  #output for eQTM panel 
  output$corTable<- renderTable( {
    
    cor_data <-
      dat_cor %>%
      filter(gene_name == input$geneInput_cor[1] & adjusted.pvalue< 0.05)%>% 
      arrange(adjusted.pvalue)%>%
      select(Gene=gene_name, 
             Probe=probe,
             Estimate=estimate,
             Pvalue=pvalue,
             Adjusted_Pvalue=adjusted.pvalue)    
  })
  
  #output for multiple probe panel 
  output$corPlot <- renderPlot({
    #From Step 3
    rawdata.for.specific.gene <-cor_test_results %>% filter(gene==input$geneInput_cor2[1])
    rawdata.for.specific.gene <-          
      rawdata.for.specific.gene[order(rawdata.for.specific.gene$adjusted.pvalue),]
    gene.expressions<-data.frame(expression=subjects_genes[,input$geneInput_cor2[1]])
    
    names.of.ordered.probes<-rawdata.for.specific.gene$probe
    ordered.probe.value<-as.data.frame(t(probes_subjects[as.character(names.of.ordered.probes),]))
    designated.dataframe.for.specific.gene<-cbind(gene.expressions,ordered.probe.value)
    number.of.picked.probes <- input$topProbes
    number.of.picked.probes <- number.of.picked.probes +1
    datacor = cor(designated.dataframe.for.specific.gene[1:number.of.picked.probes])
    corrplot(datacor, method = "color", addCoef.col="grey", type= "lower") 
    
    
    
  })  
  #output for Multiple probe panel 
  output$corPlot2 <- renderPlot({
    #From Step 3
    rawdata.for.specific.gene <-cor_test_results %>% filter(gene==input$geneInput_cor2[1])
    
    rawdata.for.specific.gene <-          
      rawdata.for.specific.gene[order(rawdata.for.specific.gene$adjusted.pvalue),]
    gene.expressions<-data.frame(expression=subjects_genes[,input$geneInput_cor2[1]])
    
    names.of.ordered.probes<-rawdata.for.specific.gene$probe
    ordered.probe.value<-as.data.frame(t(probes_subjects[as.character(names.of.ordered.probes),]))
    designated.dataframe.for.specific.gene<-cbind(gene.expressions,ordered.probe.value)
    number.of.picked.probes <- input$topProbes
    number.of.picked.probes <- number.of.picked.probes +1
    datacor = cor(designated.dataframe.for.specific.gene[1:number.of.picked.probes])
    plot_correlation(designated.dataframe.for.specific.gene,number.of.picked.probes)
    
    
  })  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
```

<!--html_preserve-->
Shiny applications not supported in static R Markdown documents

<!--/html_preserve-->
