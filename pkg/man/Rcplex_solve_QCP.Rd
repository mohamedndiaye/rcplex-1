\name{Rcplex_solve_QCP}
\alias{Rcplex_solve_QCP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Solve quadratically constrained optimization problem with CPLEX}
\description{
  Interface to CPLEX solvers for quadratically constrained linear,
  quadratic, and mixed-integer programs. The general statement of the
  problem is
  \deqn{\min \frac{1}{2}x'Qx + c'x}{min 0.5x'Qx + c'x}
  \deqn{\mathrm{s.t} Ax \leq b}{s.t Ax <= b}
  \deqn{\mathrm{and} a_i'x + x'Q_ix \leq r_i for i=1,\ldots,q}{and a_i'x + x'Q_i x <= r_i for i=1,...,q}
  \deqn{lb \leq x \leq ub}{lb <= x <= ub}

  If \code{Q==NULL} then the problem is linear, if any value of the \code{vtype}
  argument is "B" or "I" then the problem is a mixed-integer program.
  The \code{control} argument is used to set CPLEX's many parameters. See
  details. The \code{objsense} determines if the problem is a
  maximization or minimization problem. The \code{sense} argument is
  used to set the constraint directions.
}
\usage{
Rcplex_solve_QCP(cvec, Amat, bvec, Qmat = NULL, QC,
       lb = 0, ub = Inf, sense = "L", objsense = c("min", "max"), vtype
  = NULL, n = 1, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cvec}{ The linear coefficient of the objective function }
  \item{Amat}{ The constraint matrix (requires \code{ncol(Amat)==length(cvec)}) }
  \item{bvec}{ The constraints right-hand side (requires \code{length(bvec)==nrow(Amat)})}
  \item{Qmat}{ The quadratic coefficient of the objective function. If
    \code{NULL} the problem is linear. If not \code{NULL}, it must be a symmetric positive
    semidefinite matrix of size \code{length(cvec)} by \code{length(cvec)}. Default \code{NULL} }
  \item{QC}{a list with three elements: \code{QC}, \code{dir}, and
  \code{b}. The element \code{QC} is a list with the quadratic part
  \code{Q}, a matrix, and the linear part of the constraint \code{L}, a
  numeric (currently nonzero values are not supported). \code{dir} has the same meaning as argument \code{sense} and
  \code{b} as \code{bvec}.} 
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
See function \code{link[Rcplex]{Rcplex}()} for more information about
  sparse matrix representation and control arguments.
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
\references{IBM ILOG CPLEX Optimization Studio documentation}
\author{Hector Corrada Bravo and Stefan Theussl}
\seealso{\code{\link{Rcplex.close}}, \code{\link{optim}}}
\examples{
## objective function
c <- c(1, 2, 3)
Q <- matrix(c(-33, 6, 0, 6, -22, 11.5, 0, 11.5, -11), nrow = 3)

## constraints

## linear part
A <- matrix(c(-1, 1, 1, -3, 1, 1), nrow = 2)
dir <- c("L", "L")
b <- c(20, 30)

## quadratic part
QC <- list(QC = list(Q = list(diag(1, nrow = 3)), L = NULL), dir = "L", b = 1)

## bounds
ub <- c(40, Inf, Inf)

## solve
res <- Rcplex_solve_QCP(c,A, b, Q, ub = ub, QC = QC, sense = dir, objsense = "max")
print(res)

## solve MIQCP
res <- Rcplex_solve_QCP(c, A, b, Q, ub = ub, QC = QC,
                        sense = dir, objsense = "max", vtype = c("C", "I", "C"))

## quadratic and linear part
QC <- list(QC = list(Q = list(diag(1, nrow = 3)), L = list(c(3,4,-3))), dir = "L", b = 1)

## solve
res <- Rcplex_solve_QCP(c,A, b, Q, ub = ub, QC = QC, sense = dir, objsense = "max")
print(res)


Rcplex.close()
}
\keyword{optimize}

