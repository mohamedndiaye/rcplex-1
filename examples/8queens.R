## 8 Queens example from lpSolve (see the package documentation
## for more information)

library("lpSolve")
# Here's an example in which we want more than one solution to a problem
# in which all variables are binary: the 8-queens problem.
obj <- rep (1, 64)
q8  <- make.q8()
dir <- rep (c("=", "<"), c(16, 26))
rhs <- rep (1, 42)
## now calculate 4 solutions
x <- lp('max', obj, , dir, rhs, dense.const = q8, 
        all.bin = TRUE, num.bin.solns = 4)

x$solution

library("Rcplex")

dir3 <- rep (c("E", "L"), c(16, 26))
## gives a segfault IF simple_triplet_matrix not column major order
ord <- order(q8[, 2])
A <- simple_triplet_matrix(i = q8[, 1][ord], j = q8[, 2][ord],
                           v = q8[, 3][ord], ncol = length(obj), nrow = length(rhs))
y <- Rcplex(obj, A, rhs, lb = 0, ub = 1, sense = dir3, objsense = "max", vtype = "B", control = list(round = TRUE), n = NA);

