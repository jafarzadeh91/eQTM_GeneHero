Readme
================
Lisa Leung
2018-04-02

a-clust
-------

-This file includes the source code for the A-clust package, in addition, to source codes written by us to use the package in a convenient way for our application.

Step1.md and Step1.Rmd
----------------------

In this file, we explain the process of organizing data into an R-friendly format and extracting the essential part of data from the ROSMAP dataset. We then go through the multicollinearity issue in our CpG probes and how we correct it. We then assess the outliers in gene expression and methylation data. \#\# project\_data\_description.Rmd and project\_data\_description.md

-   These files provides a detailed description of the raw data set we have obtained for our project. It describes the ROSMAP data sets.

Step1\_files (Images Folder)
----------------------------

-   Step\_files will include all the images that are produced in Markdown file

Data
----

### ROSMAP Data

-   We used the ROSMAP raw data published by [synapse](https://www.synapse.org/) as the input to this step. Due to the confidentiality issues, interested scientists can directly contact synapse for downloading the data using the link below:

<https://www.synapse.org/#!Synapse:syn3219045>

### rosmap\_postprocV2.RData

-   The data sets produced here consists of the results from Step 1. The .RData folder consists of three data sets: probes\_subjects, subjects\_genes, and probes\_genes\_distance. The probes\_subjects data consists of DNA methylation probe data, subjects\_genes consists of gene expression values, and probes\_genes\_distance consists of the distances on probes and genes where 0 represents a distance over 1Mb. Since it is too large to be uploaded in github, the data set could be found in this Google Drive Folder:

<https://drive.google.com/drive/folders/1u7J2reJVtPl2IVVytpTqWqJm-76lu0OY?usp=sharing>
