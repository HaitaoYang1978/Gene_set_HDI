# read the example data
example.data<-readRDS('example.data.RDS')
X<-example.data[,-dim(example.data)[2]]
Y<-example.data[,dim(example.data)[2]]
n_x<-nrow(X)


