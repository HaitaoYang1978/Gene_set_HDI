p.adjust.wy<-function (cov, pval, N = 10000) 
{
  zz <- mvrnorm(N, rep(0, ncol(cov)), cov)
  zz2 <- scale(zz, center = FALSE, scale = sqrt(diag(cov)))
  Gz <- apply(2 * pnorm(abs(zz2), lower.tail = FALSE), 1, min)
  ecdf(Gz)(pval)
}