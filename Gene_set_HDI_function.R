Gene_set_HDI<-function(X,Y,G_Var){
  #################### Variable screening via high dimensional ordinary least squares projection (HOLP) ####################
  library(screening)
  source('HOLP.R',local = TRUE)
  x=X;y=Y
  #b=1 # b=1:test the gene 1; b=2: test the gene 2
  #if (b==1){
  #  Group_test=G_1_s0 # G_1_s0: the location of gene 1 after variable screening
  #}  
  #if (b==2){
   # Group_test=G_2_s0 # G_2_s0: the location of gene 2 after variable screening
  #}  
  G_Var=list()
  #################### statistical inference with de-sparsified LASSO estimator ####################
  
  fit.lasso <- lasso.proj(x, y,multiplecorr.method = 'WY',robust = TRUE)
  G_Var[[1]]<-which(names(fit.lasso$pval) %in% gene1_name)
  G_Var[[2]]<-which(names(fit.lasso$pval) %in% gene2_name)
  
  p_value<-list()
  P_G<-c()
  P.arta_Cauchy<-c()
  P.arta_Omnibus<-c()
  Results<-list()
  for (i in 1:length(G_Var)) {
    p_value[[i]]=fit.lasso$pval[G_Var[[i]]] # the p-values of all SNPs within the gene to be tested.
    bhat=fit.lasso$bhat[G_Var[[i]]]
    
    ########### Inference under the DSSM model assumption by minimum p-value approach ###########
    
    cov_G=fit.lasso$beta.cov[G_Var[[i]],G_Var[[i]]]
    P_G[i]=min(p.adjust.wy(cov=cov_G, pval=p_value[[i]]))
    
    #################### Inference under the CWSM model assumption by iART-A ####################
    
    source('DOT.ART.P.R',local = TRUE) #Decorrelation by orthogonal transformation (DOT)
    source('ART.A.R',local = TRUE) #Adaptative augmented rank truncation (ART-A)
    
    #################### Improved adaptative augmented rank truncation (iART-A) ####################
    L<-length(p_value[[i]])
    k1 <- 2
    k2 <-L
    P.arta<-c()
    for ( k in k1:k2){
      P.arta[which(k== k1:k2)]<- ART.A(P, k, L)[1]
    }
    
    P.arta_Cauchy[i]<-ACAT(P.arta)
    d<-L-1
    P.arta_Cauchy[i]<-ifelse(P.arta_Cauchy[i]==1,P.arta_Cauchy[i]<-1-1/d,P.arta_Cauchy[i])#replace 1 by 1-1/d, where d is the number of p-values combined by ACAT
    
    #################### Omnibus test based on MinP and iART-A ####################
    
    P.arta_Omnibus[i]<-ACAT(c(P_G,P.arta_Cauchy))
    
    #################### Output results into a exceldata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg== file ####################
    
  }
  names(p_value)<-c('gene1','gene2')
  names(P_G)<-c('gene1','gene2')
  names(P.arta_Cauchy)<-c('gene1','gene2')
  names(P.arta_Omnibus)<-c('gene1','gene2')

  Results<-list(p_value,P_G,P.arta_Cauchy,P.arta_Omnibus)
  names(Results)[1]<-c('p_value')
  names(Results)[2]<-c('MinP')
  names(Results)[3]<-c('iART_A')
  names(Results)[4]<-c('Min_O') 
 
  return(Results)
  
}