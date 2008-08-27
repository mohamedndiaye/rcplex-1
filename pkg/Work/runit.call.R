test.lp1 <- function()
{
  # A linear program (this is lpex1.c in the CPLEX examples)
cvec <- c(1,2,3);
Amat <- matrix(c(-1,1,1,-1,3,-1),byrow=TRUE,nc=3);
bvec <- c(20,-30);
ub <- c(40,Inf,Inf);

  res <- Rcplex(cvec,Amat,bvec,ub=ub,objsense="max",sense=c('L','G'));
  checkEquals(res$xopt,c(40,17.5,42.5));
  checkEquals(res$obj,202.5);
  checkEquals(res$status,1);
}

test.lp2 <- function()
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
x0 <- runif(n);
b <- as.matrix(Amat) %*% x0;
cvec <- rnorm(n);

res <- Rcplex(cvec,Amat,b,sense='E',control=list(method=4));
print(res$obj);print(crossprod(cvec,x0));
checkEquals(res$obj,as.numeric(crossprod(cvec,x0)),check.attributes=FALSE);
checkEquals(res$status,1);
}

test.qp1 <- function()
  {
# A quadratic problem (this is qpex1.c in the CPLEX examples)
cvec <- c(1,2,3);
Qmat <- matrix(c(-33,6,0,
                  6,-22,11.5,
                  0,11.5,-11),
                byrow=TRUE,
                nc=3);
Amat <- matrix(c(-1,1,1,
                  1,-3,1),
               byrow=TRUE,nc=3);
bvec <- c(20,30);
ub <- c(40,Inf,Inf);

res <- Rcplex(cvec,Amat,bvec,Qmat,ub=ub,objsense="max"); print(res);
checkEquals(res$xopt,c(0.1391,0.5985,0.8984),check.attr=FALSE,tol=1e-4);
checkEquals(res$obj,2.0156,tol=1e-4);
checkEquals(res$status,1);    
}

test.milp1 <- function()
  {
# A mixed integer linear program (mipex1.c in the CPLEX examples)
cvec <- c(1,2,3,1);
Amat <- matrix(c(-1,1,1,10,
                  1,-3,1,0,
                  0,1,0,-3.5),
               byrow=TRUE, nc=4);
bvec <- c(20,30,0);
lb <- c(0,0,0,2);
ub <- c(40,Inf,Inf,3);
vtype <- c(rep("C",3),"I");

res <- Rcplex(cvec,Amat,bvec,lb=lb,ub=ub,sense=c("L","L","E"),objsense="max",vtype=vtype);
checkEquals(res$xopt,c(40,10.5,19.5,3));
checkEquals(res$obj,122.5);
checkEquals(res$status,101);
}

test.miqp1 <- function()
  {
# A mixed integer quadratic program
cvec <- c(1,2,3,1);
Qmat <- matrix(c(-33,6,0,0,
                  6,-22,11.5,0,
                  0,11.5,-11,0,
                  0,0,0,0),
               byrow=TRUE, nc=4);
Amat <- matrix(c(-1,1,1,10,
                  1,-3,1,0,
                  0,1,0,-3.5),
               byrow=TRUE, nc=4);
bvec <- c(20,30,0);
ub <- c(40,Inf,Inf,3);
vtype <- c(rep("C",3),"I");

res <- Rcplex(cvec,Amat,bvec,Qmat=Qmat,ub=ub,sense=c("L","L","E"),objsense="max",vtype=vtype);
checkEquals(res$xopt,c(0.0303,0,0.2727,0),tol=1e-3);
checkEquals(res$obj,0.4242,tol=1e-3);
checkEquals(res$status,101);
}

