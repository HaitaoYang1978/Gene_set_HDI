library(MASS)

#generate multivariate normal covariance matrix sigma for b1
G1<-readRDS("G1.RDS")
G2<-readRDS("G2.RDS")
G_N_1<-dim(G1)[2]
G_N_2<-dim(G2)[2]
X<-cbind(G1,G2)
set.seed(5*gg+1)
n_id<-sample(1:nrow(X),n_x)
set.seed(5*gg)
epsilon=rnorm(n_x,0,1)
X_G1=X[n_id,1:G_N_1]
X_G2=X[n_id,(G_N_1+1):ncol(X)]
b1=as.matrix(c(runif(n_c_1,c_beta_1,c_beta_1),rep(0,G_N_1-n_c_1)))
set.seed(5*gg)
b1<-sample(b1)
b2=as.matrix(c(runif(n_c_2,c_beta_2,c_beta_2),rep(0,G_N_2-n_c_2)))
set.seed(5*gg)
b2<-sample(b2)
X_G1<-as.matrix(X_G1)
X_G2<-as.matrix(X_G2)
X=cbind(X_G1,X_G2)
X=scale(X,center = TRUE,scale = TRUE)



G_1_Test=c(1:G_N_1)
G_2_Test=c((G_N_1+1):(G_N_1+G_N_2))
#Size or Power
Group_test=G_1_Test

