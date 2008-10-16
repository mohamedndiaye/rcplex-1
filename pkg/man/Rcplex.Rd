\name{Rcplex}
\alias{Rcplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Solve optimization problem with CPLEX}
\description{
  Interface to CPLEX solvers for linear quadratic and (linear or
  quadratic) mixed-integer programs. The general statement of the
  problem is
  \deqn{\min \frac{1}{2}x'Qx + c'x}{min 0.5x'Qx + c'x}
  \deqn{\mathrm{s.t} Ax \leq b}{s.t Ax <= b}
  \deqn{lb \leq x \leq ub}{lb <= x <= ub}

  If \code{Q==NULL} then the problem is linear, if any value of the \code{vtype}
  argument is "B" or "I" then the problem is a mixed-integer program.
  The \code{control} argument is used to set CPLEX's many parameters. See
  details. The \code{objsense} determines if the problem is a
  maximization or minimization problem. The \code{sense} argument is
  used to set the constraint directions.
}
\usage{
Rcplex(cvec, Amat, bvec, Qmat = NULL,
       lb = 0, ub = Inf, control = list(),
       objsense = c("min", "max"), sense = "L", vtype = NULL, n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cvec}{ The linear coefficient of the objective function }
  \item{Amat}{ The constraint matrix (requires \code{ncol(Amat)==length(cvec)}) }
  \item{bvec}{ The constraints right-hand side (requires \code{length(bvec)==nrow(Amat)})}
  \item{Qmat}{ The quadratic coefficient of the objective function. If
    \code{NULL} the problem is linear. If not \code{NULL}, it must be a symmetric positive
    semidefinite matrix of size \code{length(cvec)} by \code{length(cvec)}. Default \code{NULL} }
  \item{lb}{ Lower bound on the problem variables. If
    \code{length(lb)==1} then \code{lb} is the lower bound of all
    variables. Otherwise, \code{length(lb)==length(cvec)}. Set
    \code{lb=-Inf} to have no lower bound. Default 0.}
  \item{ub}{ Upper bound on the problem variables. See \code{lb} for
    further details. Default \code{Inf}.}
  \item{control}{ A list of CPLEX parameters. See *Details* }
  \item{objsense}{ Either \code{"max"} or \code{"min"}, determines the optimization
    direction. Default \code{"min"}}
  \item{sense}{ The direction of the inequality in each
    constraint. If \code{length(sense)==1} then the same value is taken
    for each constraint. Can be one of \code{"L"} (less than or equal),
    \code{"G"} (reater than or equal) or \code{"E"} (equal). Requires
    \code{length(sense)==length(bvec)}. Default \code{"L"}. }
  \item{vtype}{ Determines the type of each problem variable. Can be one
    of \code{"C"} (continuous), \code{"I"} (integer) or \code{"B"} (binary). If
    \code{length(vtype)==1} the same value is taken for all
    variables. Otherwise, requires
    \code{length(vtype)==length(ctype)}. Default \code{"C"}.}
  \item{n}{ Determines the maximal number of solutions the solver should
    return in case of an MIP with more than one solution at
    optimum. If CPLEX should search for "all" solutions then
    \code{n} has to be set to \code{NA}. In CPLEX this is also called
    populating the solution pool. The parameters \code{solnpoolagap},
    \code{solnpoolgap}, and \code{solnpoolintensity} influence the
    search for multiple solutions (see also the \code{control}
    argument below for details). Available from CPLEX 11.0 on. \code{Rcplex()}
    raises a warning if an older version of CPLEX is used and \code{n>1}. Default \code{1}.} 
}
\details{
  Matrices \code{A} and \code{C} may be sparse matrices from a class in
  the hierarchy defined by the \pkg{Matrix} package. In that case, the
  internal casting functions are used to create the proper data
  structures to pass to CPLEX, which is similar to the column-major
  storage mode defined by the \code{dgCMatrix-class} defined by the
  \pkg{Matrix} package.

  We also provide a simple S3-style class for sparse matrices
  \code{simple_triplet_matrix}, as used in the \pkg{relations}
  package. Matrices \code{A} and \code{C} can be objects of this class.
  See the \code{examples} for example
  usage. \code{simple_triplet_matrix} objects MUST be in column-major
  order.
  
  The \code{control} argument can be used to set CPLEX's many
  parameters, including the particular algorithm used for solving the given
  problem. See the \cite{ILOG CPLEX Parameters} guide for further
  details. The following parameters are supported:

  \describe{
    \item{trace:}{Turn CPLEX output on (1) or off(0). Default 1.}
    \item{maxcalls:}{Number of calls to the CPLEX optimizer before
      license is released. Set to 1 to get a new license on every call
      to Rcplex. Can be any positive number. Default 500.}
    \item{method:}{Algorithm to use (Default 0):
      \describe{
	\item{0:}{Automatic: CPLEX chooses algorithm automatically}
	\item{1:}{Primal Simplex}
	\item{2:}{Dual Simplex}
	\item{3:}{Network Simplex}
	\item{4:}{Barrier}
      }
    }
    \item{preind:}{Turn presolver on (1) or off (0). Default 1.}
    \item{aggind:}{Limit on the number of applications of the
      aggregator. Possible Values:
      -1 (automatic), 0 (do not use), any positive integer}
    \item{itlim:}{Maximum number of simplex iterations. Can be any
      nonnegative number. Default 1e8.}
    \item{epagap:}{Absolute MIP optimality gap tolerance. Can be any
      nonnegative number. Default 1e-6.}
    \item{epgap:}{Relative MIP optimality gap tolerance. Can be any
      nonnegative number. Default 1e-4.}
    \item{tilim:}{Time limit in seconds of call to optimizer. Can be any
      nonnegative number. Default 1e75.}
    \item{disjcuts:}{Indicator for disjunctive cuts used in MIP solver.
      Must be in -1:3. Default 0 (automatic).}
    \item{mipemphasis:}{Indicator for MIP solver emphasis. Must be in
      0:4. Default 0 (balance optimality and feasibility)}
    \item{cliques:}{Indicator for clique cuts in MIP solver. Must be in
      -1:2. Default 0 (automatic)}
    \item{nodesel:}{Node selection strategy in MIP solver. Must be in
      0:3. Default 1 (best-bound search).}
    \item{probe:}{Probe level in MPI solver. Must be -1:3. Default 0
      (automatic)}
    \item{varsel:}{Variable selection strategy in MIP solver. Must be in
      -1:4. Default 0 (choose best method automatically).}
    \item{flowcovers:}{Indicator for flowcover cuts in MIP solver. Must
      be in -1:2. Default 0 (automatic).}
    \item{solnpoolagap:}{Sets an absolute tolerance on the objective
      value for the solutions in the solution pool. Can be any
      nonnegative real number. Ignored in versions < 11.0 of CPLEX. Default 0}
    \item{solnpoolgap:}{Sets a relative tolerance on the objective
      value for the solutions in the solution pool. Can be any
      nonnegative real number. Ignored in versions < 11.0 of CPLEX. Default 0}
    \item{solnpoolintensity:}{Controls the trade-off between the number
      of solutions generated for the solution pool and the amount of
      time and memory consumed. Must be in 0:4. Ignored in versions <
      11.0 of CPLEX. Default 0 (automatic).}
    \item{round:}{Flag indicating if integer solutions for MIPs should be
      rounded before returning. In some cases, CPLEX returns slightly
      infeasible integer solutions. Setting this option to 1 ensures
      that the returned solution is integral by rounding. Default 0 (no rounding).}
  }
}
\value{
  Returns a list with the following components, or, if \code{n > 1} a list
  of length equal to the number of optimal solutions containing the
  following components for each solution:
  \item{xopt}{Values of problem variables at optimum.}
  \item{obj }{Value of objective function at optimum.}
  \item{status}{Solution status. See
    CPLEX documentation for meaning of status codes.}
  \item{extra}{List with extra information about solution with components
    \describe{
      \item{slack:}{Values of slack variables for inequality
	constraints.}
      \item{nodecnt:}{(IF MIP PROBLEM) Number of nodes in the search tree evaluated}
      \item{lambda:}{(IF NOT MIP PROBLEM) Values of dual variables at optimum}
    }
  }
}
\references{ILOG CPLEX User's Guide }
\author{Hector Corrada Bravo and Stefan Theussl }
\seealso{\code{\link{Rcplex.close}}, \code{\link{optim}}}
\examples{
## A linear program (this is lpex1.c in the CPLEX examples)
cvec <- c(1,2,3)
Amat <- matrix(c(-1,1,1,-1,3,-1),byrow=TRUE,nc=3)
bvec <- c(20,-30)
ub <- c(40,Inf,Inf)

res <- Rcplex(cvec,Amat,bvec,ub=ub,objsense="max",sense=c('L','G'))
print(res)

## A linear program with random data
## use the barrier method
n = 20; m = 25
nnz <- trunc(.2 * m * n)

## entries in simple_triplet_matrix clas
##  *must* be in column major order
nnz <- sort(sample(m*n,nnz,replace=FALSE)-1)
Amat <- simple_triplet_matrix(
             i = (nnz \%\% m) + 1,
             j = trunc(nnz/m) + 1,
             v = rnorm(nnz),
             nrow=m,ncol=n)

x0 <- runif(n)
b <- as.matrix(Amat) \%*\% x0
cvec <- rnorm(n)

res <- Rcplex(cvec,Amat,b,sense='E',control=list(method=4))
print(res)

## A quadratic problem (this is qpex1.c in the CPLEX examples)
cvec <- c(1,2,3)
Qmat <- matrix(c(-33,6,0,
                  6,-22,11.5,
                  0,11.5,-11),
                byrow=TRUE,
                nc=3)
Amat <- matrix(c(-1,1,1,
                  1,-3,1),
               byrow=TRUE,nc=3)
bvec <- c(20,30)
ub <- c(40,Inf,Inf)

res <- Rcplex(cvec,Amat,bvec,Qmat,ub=ub,objsense="max")
print(res)

## A mixed integer linear program (mipex1.c in the CPLEX examples)
cvec <- c(1,2,3,1)
Amat <- matrix(c(-1,1,1,10,
                  1,-3,1,0,
                  0,1,0,-3.5),
               byrow=TRUE, nc=4)
bvec <- c(20,30,0)
lb <- c(0,0,0,2)
ub <- c(40,Inf,Inf,3)
vtype <- c(rep("C",3),"I")

res <- Rcplex(cvec,Amat,bvec,lb=lb,ub=ub,sense=c("L","L","E"),
              objsense="max",vtype=vtype)
print(res)

## A mixed integer quadratic program
cvec <- c(1,2,3,1)
Qmat <- matrix(c(-33,6,0,0,
                  6,-22,11.5,0,
                  0,11.5,-11,0,
                  0,0,0,0),
               byrow=TRUE, nc=4)
Amat <- matrix(c(-1,1,1,10,
                  1,-3,1,0,
                  0,1,0,-3.5),
               byrow=TRUE, nc=4)
bvec <- c(20,30,0)
ub <- c(40,Inf,Inf,3)
vtype <- c(rep("C",3),"I")

res <- Rcplex(cvec,Amat,bvec,Qmat=Qmat,ub=ub,sense=c("L","L","E"),
              objsense="max",vtype=vtype)
print(res)
Rcplex.close()
}
\keyword{optimize}

