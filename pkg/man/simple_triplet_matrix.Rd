\name{simple_triplet_matrix}
\alias{simple_triplet_matrix}
\title{Create sparse matrix with triplet representation}
\description{
  Create a sparse matrix, with entries stored by triplets (row,col,entry).
}
\usage{
simple_triplet_matrix(i, j, v, nrow = max(i), ncol = max(j))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{Row indices of entries}
  \item{j}{Column indices of entries}
  \item{v}{Entries}
  \item{nrow}{Number of rows in matrix}
  \item{ncol}{Number of columns in matrix}
}
\details{
  For usage in \code{Rcplex}, entries must be given in column major
  order, i.e. component \code{j} must be in ascending order
}
\value{
  Returns \code{S3}-style object of class \code{simple_triplet_matrix}
  \item{i}{Row indices of entries}
  \item{j}{Column indices of entries}
  \item{v}{Entries}
  \item{nrow}{Number of rows in matrix}
  \item{ncol}{Number of columns in matrix}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Hector Corrada Bravo <hcorrada@gmail.com>, uses code from the
  \code{relations} package of Kurt Hornik and David Meyer}
\examples{
## Define matrix
## 1 0 3
## 0 1 2
## 3 2 0
mat <- simple_triplet_matrix(
   i=c(1,3,2,3,1,3),
   j=c(1,1,2,2,3,3),
   v=c(1,3,1,2,3,2));
print(mat);

# coerce to matrix
as.matrix(mat);

# transpose
as.matrix(t(mat));

# indexing
mat[1,1];
mat[1,];
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
