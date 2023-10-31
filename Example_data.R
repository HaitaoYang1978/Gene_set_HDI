# read the example data
example.data<-readRDS('example.data.RDS')
X<-example.data[,-dim(example.data)[2]]
Y<-example.data[,dim(example.data)[2]]
X=scale(X,center = TRUE,scale = TRUE)
G_N_1<-168 # The number of SNPs in Gene 1
G_N_2<-823 # The number of SNPs in Gene 2
n_x<-nrow(X)
X_G1=X[,1:G_N_1]
X_G2=X[,(G_N_1+1):ncol(X)]

X_G1<-as.matrix(X_G1)
X_G2<-as.matrix(X_G2)
G_1_Test=c(1:G_N_1)
G_2_Test=c((G_N_1+1):(G_N_1+G_N_2))

