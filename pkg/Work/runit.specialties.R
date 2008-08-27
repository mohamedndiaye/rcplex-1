test.round <- function()
  {
    # A linear program with random data
# use the barrier method
n = 20; m = 25;
nnz <- trunc(.2 * m * n);
nnz <- sort(sample(m*n,nnz,replace=FALSE)-1);
Amat <- simple_triplet_matrix(
             i = (nnz %% m) + 1,
             j = trunc(nnz/m)+1,
             v = rnorm(nnz),
             nrow=m,ncol=n);

x0 <- 100*runif(n);
intvars <- 1:(n/2)
x0[intvars] <- round(x0[intvars])
  
b <- as.matrix(Amat) %*% x0;
cvec <- rnorm(n);
cvec[intvars] <- 0

vtype <- rep("C",n)
vtype[intvars] <- "I"
  
res <- Rcplex(cvec,Amat,b,sense='E',vtype=vtype,control=list(preind=0,round=1));
print(res$obj);print(crossprod(cvec,x0));
checkEquals(res$obj,as.numeric(crossprod(cvec,x0)),check.attributes=FALSE);
checkEquals(res$status,101);
checkEquals(round(res$xopt[intvars]),res$xopt[intvars])
  }
