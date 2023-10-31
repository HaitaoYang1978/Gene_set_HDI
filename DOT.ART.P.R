## the number of statistics

Z <- as.vector(fit.lasso$beta.statistic[Group_test])
L <- length(Z)
## transform Z-scores to chi-squares
y <- Z^2
## calculate TQ-statistic
yy <- y %*% rep(1, L)

S<-fit.lasso$beta.cor[Group_test,Group_test]
ee <- eigen(S); eivec <- ee$vectors; eigva <- ee$values
pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
## calculate decorreated statistics
x <- (Z %*% pc)^2
px <- 1-pchisq(x, df=1)
P <- sort(px)
