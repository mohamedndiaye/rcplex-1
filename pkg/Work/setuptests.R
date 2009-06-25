pkg.dir <- file.path("..");
src.files <- list.files(file.path(pkg.dir,"R"),pattern=".*\\.R$")
invisible(lapply(src.files,function(filenm) source(file.path(pkg.dir,"R",filenm))))
dyn.load(file.path(pkg.dir,"src",paste("Rcplex",.Platform$dynlib.ext,sep="")))

require("slam")
cvec <- c(1,2,3);
Amat <- matrix(c(-1,1,1,-1,3,-1),byrow=TRUE,nc=3);
bvec <- c(20,-30);
ub <- c(40,Inf,Inf);
Qmat <- matrix(c(-33,6,0,
                  6,-22,11.5,
                  0,11.5,-11),
                byrow=TRUE,
                nc=3);

stmAmat <- structure(list(i=c(1,2,1,2,1,2),
                          j=c(1,1,2,2,3,3),
                          v=c(-1,-1,1,3,1,-1),
                          nrow=2,ncol=3),
                     class="simple_triplet_matrix");

