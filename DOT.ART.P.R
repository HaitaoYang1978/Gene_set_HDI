## the number of statistics

Z <- as.vector(fit.lasso$beta.statistic[G_Var[[i]]])
L <- length(Z)
## transform Z-scores to chi-squares
y <- Z^2
## calculate TQ-statistic
yy <- y %*% rep(1, L)

S<-fit.lasso$beta.cor[G_Var[[i]],G_Var[[i]]]
ee <- eigen(S); eivec <- ee$vectors; eigva <- ee$values
pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
## calculate decorreated statistics
x <- (Z %*% pc)^2
px <- 1-pchisq(x, df=1)
P <- sort(px)
