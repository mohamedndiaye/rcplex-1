## Example from cplex documentation:
## Maximize
## obj: x1 + 2 x2 + 3 x3 + [ - 33 x1 ^2 + 12 x1 * x2 - 22 x2 ^2 + 23 x2 * x3
##      - 11 x3 ^2 ] / 2
## Subject To
##  c1: - x1 + x2 + x3 <= 20
##  c2: x1 - 3 x2 + x3 <= 30
##  q1: [ x1 ^2 + x2 ^2 + x3 ^2 ] <= 1
## Bounds
##  0 <= x1 <= 40
## End

library("Rcplex")

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

Rcplex:::Rcplex_solve_LP(c,A, b, ub = ub, sense = dir, objsense = "max")

Rcplex:::Rcplex_solve_QP(c,A, b, Q, ub = ub, sense = dir, objsense = "max")

Rcplex:::Rcplex_solve_QCP(c,A, b, ub = ub, QC = QC, sense = dir, objsense = "max")

Rcplex:::Rcplex_solve_QCP(c,A, b, Q, ub = ub, QC = QC, sense = dir, objsense = "max")


## quadratic and linear part
QC <- list(QC = list(Q = list(diag(1, nrow = 3)), L = list(c(3,4,-3))), dir = "L", b = 1)

Rcplex:::Rcplex_solve_QCP(c,A, b, Q, ub = ub, QC = QC, sense = dir, objsense = "max")


q(save = "no")

