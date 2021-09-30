# CIS4930 Project Repository Fall 2021
## Team Members

Name | Email | Github
---- | ----- | ------
Michael Mezzina | michaelmez39@gmail.com | github.com/michaelmez39
Zhiqing Qu | zhiqingqu@ufl.edu | github.com/ZhiqingQu

## Research Question
Identify genetic differences between lung organoids infected with covid19 and those not infected.

## Dataset
[Covid 19 Datasets](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL29320)<br>
[Our Current Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162323)<br>

## Setup
1. Download R, RStudio
3. Download files "GSE162323_slam_inf_params.txt.gz" and "GSE162323_slam_inf_params.txt.gz" from [the dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162323)
4. Unzip these files.
5. Make folder "data" and move files the "data" folder
6. Open this project in R Studio
7. Install "dplyr" and "data.table" with cran e.g. `R> packages.install("dplyr")`
8. Install "EnhancedVolcano", "DESeq2", and "apeglm" with bioconductor. `R>BiocManager::install("DESeq2", update= FALSE)`
