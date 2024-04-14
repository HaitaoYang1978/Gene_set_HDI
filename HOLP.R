library(screening)
set.seed(gg)
train_id<-sample(1:n_x, n_x*0.5)
d=floor(dim(X[train_id,])[1]/log(dim(X[train_id,])[1]))
model11=screening(X[train_id,], Y[train_id], method = "holp", num.select =d,family = "gaussian", ebic = FALSE)
s0<-sort(model11$screen)
Xt<-X[-train_id,]
Yt<-Y[-train_id]
X<-Xt[,s0]
Y<-Yt










