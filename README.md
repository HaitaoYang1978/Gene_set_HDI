# Gene_set_HDI
## Description
Our proposed procedure aims to detect the set-based genetic variation (e.g., SNP-set, gene-set) associated with complex disease trait under a high-dimensional inference framework. Furthermore, we propose an omnibus testing procedure that employs a robust and powerful p-value combination method to enhance the power of set-based genetic variation association.
## Dependency
Our proposed procedure relies on the package `glmnet`  `future.apply` `MASS` `mvtnorm` `devtools` `ACAT` `screening`. All required package can be installed automatically in our procedure. In addition, all corresponding R functions of de-sparsified lasso method are in the *'./de-sparsified lasso functions'* folder.
## Usage
Gene_set_HDI(X,Y,G_Var=G_var)
## Argument
X        Predictors. An n by p matrix, p is the number of predictors. The predictors can be discrete (e.g., SNP 
         genotype) or continuous (e.g., RNA-seq).
Y        A response vector that can be discrete or continuous.
G_Var    The location of genes or pathway
## Value
p_value   The original p values of SNPs within some a gene. A list. 
MinP      The p value of each gene via the maximum statistics distribution. A vector.
iART_A     The p value of each gene via the improved adaptative augmented rank truncation (iART-A). A vector.
Min_O      The p value of each gene via Omnibus test based on MinP and iART-A
## Author(s)
Haitao Yang and Fuzhao Chen. We appreciate your feedback under issues of this repository.
## References
1. Fan, Jianqing, and Jinchi Lv. "Sure independence screening for ultrahigh dimensional feature space." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.5 (2008): 849-911.
2. Wang X, Leng C. High dimensional ordinary least squares projection for screening variables, Journal of the Royal Statistical Society: Series B (Statistical Methodology) 2016;78:589-611.3. Li, Gaorong, et al. "Robust rank correlation based screening." The Annals of Statistics 40.3 (2012): 1846-1877.
3. Liu Y, Xie J. Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures, Journal of the American Statistical Association 2019:1-18.
4. Vsevolozhskaya OA, Hu F, Zaykin DV. Detecting weak signals by combining small P-values in genetic association studies, Frontiers in Genetics 2019;10:1051.
5. Zhang CH, Zhang SS. Confidence intervals for low dimensional parameters in high dimensional linear models, Journal of the Royal Statistical Society: Series B (Statistical Methodology) 2014;76:217-242.
6. Van de Geer S, Bühlmann P, Ritov Ya et al. On asymptotically optimal confidence regions and tests for high-dimensional models, The Annals of Statistics 2014;42:1166-1202.

## Tutorial based on a toy data of SNP genotype
A simple working toy data, a 600×992 matrix consisting of two genes and a phenotype. There are 168 SNPs in gene 1 and 823 SNPs in gene 2. the gene2 is set to be associated with the phenotype dependent on the cumulative weak signal model (CWSM).
### Set Working Directory and random seed
```
#setwd ('~./Gene_set_HDI/')
gg=2 # set the random seed.
```
### load all the required packages

```
if(!require(glmnet)){ 
  install.packages("glmnet")
  library(glmnet)
  }
if(!require(future.apply)){ 
  install.packages("future.apply") 
  library(future.apply)
  }
if(!require(MASS)){ 
  install.packages("MASS")
  library(MASS)
  }
if(!require(mvtnorm)){ 
  install.packages("mvtnorm") 
  library(mvtnorm)
  }
if(!require(devtools)){ 
  install.packages("devtools")
  library(devtools) }
if(!require(ACAT)){ 
  library(devtools)
  devtools::install_github("yaowuliu/ACAT") 
  library(ACAT)
}
if(!require(screening)){ 
  library(devtools)
  devtools::install_github('wwrechard/screening') 
  library(screening)
}
```
### Source all de-sparsified lasso functions
```
files.sources = list.files('./de-sparsified lasso functions')
sapply(paste0("./de-sparsified lasso functions/",files.sources), source)
```
### Load the example data
`source('Example_data.R')`
### Mapping the SNPs to the genes
```
# The real SNPs data can be mapped to genes by using biomaRt package 
total_SNP_name<-colnames(X)
gene1_name<-colnames(X)[1:168]
gene2_name<-colnames(X)[169:991]
G1_location<-which(total_SNP_name %in% gene1_name)
G2_location<-which(total_SNP_name %in% gene2_name)
G_Var<-list(G1_location,G2_location)
names(G_Var)<-c('gene1','gene2')
```
### Testing the genes
```
source('Gene_set_HDI_function.R')
Gene_set_HDI(X,Y,G_Var=G_var)
```
```
$p_value
$p_value$gene1
  rs970973_T rs10489142_G  rs1029322_T  rs1750837_G rs11805932_A rs10864304_G 
   0.6139564    0.8166387    0.9421041    0.2265243    0.5760912    0.7098036 
 rs7529399_C rs11120973_A  rs4908656_C rs11587479_C   rs228729_A 
   0.8630459    0.6999327    0.2145637    0.8109494    0.3231438 

$p_value$gene2
 rs7818488_A  rs7002853_C rs12545633_G rs10104274_A  rs7843496_T rs17065935_G 
 0.953658045  0.093367059  0.873342578  0.850507019  0.758513213  0.015570106 
rs13279490_T  rs7828513_T  rs2406985_A  rs4523281_T rs10503214_C  rs2720797_T 
 0.022909902  0.024747093  0.825634287  0.262280471  0.091181573  0.190369527 
 rs7843888_A  rs2623755_A rs12550117_C rs17744884_C  rs2204223_G rs17068473_T 
 0.203814360  0.535521641  0.101529279  0.317280113  0.229797994  0.162634523 
 rs2688388_C rs17068727_T  rs2554704_C  rs7018112_T rs17069257_C  rs4875308_A 
 0.597307763  0.378462675  0.117424826  0.863390180  0.099102129  0.929097651 
 rs1714717_G  rs4404937_A  rs2407940_C  rs1504757_T  rs7845852_A  rs7002661_G 
 0.525100935  0.965343674  0.715866145  0.223048746  0.879722120  0.014367339 
rs10107472_A rs11985125_T rs17069985_G  rs4875101_A rs11136742_T  rs2617086_A 
 0.233928318  0.013282781  0.843608770  0.704994131  0.007907595  0.199173114 
rs11781113_A   rs893465_T  rs3886808_G  rs6558954_T  rs4875441_A 
 0.004365529  0.959935885  0.169825641  0.753136567  0.387376981 
$MinP
 gene1  gene2 
0.9310 0.1611 
$iART_A
      gene1       gene2 
0.826674959 0.002654436 
$Min_O
     gene1      gene2 
0.90055371 0.01101687
```
When the significant level is 0.05 , The p value (Min_O) of gene2 via Omnibus test based on MinP and iART-Athe gene2 is significant, which is consistent with the setup of example data. However, The p value (MinP) of gene2 via the maximum statistics distribution is not significant.


