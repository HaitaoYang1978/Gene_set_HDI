#setwd('~./Gene_set_HDI/')
gg=2 # set the random seed.

#################### load all the required packages ####################
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

#################### Source all de-sparsified lasso functions ####################

files.sources = list.files('./de-sparsified lasso functions')
sapply(paste0("./de-sparsified lasso functions/",files.sources), source)

#################### Load the example data ####################

source('Example_data.R')

#################### Variable screening via high dimensional ordinary least squares projection (HOLP) ####################
library(screening)
source('HOLP.R')
x=X;y=Y
b=2 # b=1:test the gene 1; b=2: test the gene 2
if (b==1){
Group_test=G_1_Test # G_1_Test: the location of gene 1 after variable screening
}  
if (b==2){
Group_test=G_2_Test # G_2_Test: the location of gene 2 after variable screening
}  
G_Var=Group_test
#################### statistical inference with de-sparsified LASSO estimator ####################

fit.lasso <- lasso.proj(x, y,multiplecorr.method = 'WY',robust = TRUE)
p_value=fit.lasso$pval[Group_test] # the p-values of all SNPs within the gene to be tested.
bhat=fit.lasso$bhat[Group_test]

########### Inference under the DSSM model assumption by minimum p-value approach ###########

cov_G=fit.lasso$beta.cov[Group_test,Group_test]
P_G=min(p.adjust.wy(cov=cov_G, pval=p_value))

#################### Inference under the CWSM model assumption by iART-A ####################

source('DOT.ART.P.R') #Decorrelation by orthogonal transformation (DOT)
source('ART.A.R') #Adaptative augmented rank truncation (ART-A)

#################### Improved adaptative augmented rank truncation (iART-A) ####################
L<-length(p_value)
k1 <- 2
k2 <-L
P.arta<-c()
for ( k in k1:k2){
P.arta[which(k== k1:k2)]<- ART.A(P, k, L)[1]
}

P.arta_Cauchy<-ACAT(P.arta)
d<-L-1
P.arta_Cauchy<-ifelse(P.arta_Cauchy==1,P.arta_Cauchy<-1-1/d,P.arta_Cauchy)#replace 1 by 1-1/d, where d is the number of p-values combined by ACAT

#################### Omnibus test based on MinP and iART-A ####################

P.arta_Omnibus<-ACAT(c(P_G,P.arta_Cauchy))
 
#################### Output results into a excel file ####################
Results<-cbind(P_G,P.arta_Cauchy,P.arta_Omnibus)
Results
