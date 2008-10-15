
test.multsols <- function(){

  ## borrowed from lpSolve to construct constraint matrix for 8queens problem
  q8foo <- function () {
    chess <- matrix(1:64, 8, 8, byrow = T)
    row.const <- cbind(rep(1:8, each = 8), c(t(chess)), 1)
    col.const <- cbind(rep(9:16, each = 8), c(chess), 1)
    const.ctr <- 17
    chess <- matrix(1:64, 8, 8, byrow = T)
    row.chess <- row(chess)
    col.chess <- col(chess)
    d1.const <- NULL
    rplusc <- row.chess + col.chess
    for (i in 3:15) {
      d1.const <- rbind(d1.const, cbind(const.ctr, chess[rplusc == 
                                                         i], 1))
      const.ctr <- const.ctr + 1
    }
    start <- seq(49, 1, by = -8)
    for (i in start) {
      d1.const <- rbind(d1.const, cbind(const.ctr, seq(i, 64, 
                                                       by = 9), 1))
      const.ctr <- const.ctr + 1
    }
    for (i in 2:7) {
      d1.const <- rbind(d1.const, cbind(const.ctr, seq(i, chess[9 - 
                                                                i, 8], by = 9), 1))
      const.ctr <- const.ctr + 1
    }
    rbind(row.const, col.const, d1.const)
  }
  
  obj <- rep (1, 64)
  q8 <- q8foo()
  dir <- rep (c("E", "L"), c(16, 26))
  rhs <- rep (1, 42)
  
  ## column major order for simple_triplet matrix
  ord <- order(q8[, 2])
  constr <- simple_triplet_matrix(i = q8[, 1][ord], j = q8[, 2][ord],
                                  v = q8[, 3][ord], ncol = length(obj),
                                  nrow = length(rhs))

  ## solve problem max 1 solution 
  n_one <- Rcplex(obj, constr, rhs, lb = 0, ub = 1, sense = dir,
                  objsense = "max", vtype = "B",
                  control = list(round = TRUE), n = 1L)


  ## solve the roblem return max 4 solutions
  n_four <- Rcplex(obj, constr, rhs, lb = 0, ub = 1, sense = dir,
                   objsense = "max", vtype = "B",
                   control = list(round = TRUE), n = 4L)
  
  ## find all solutions
  n_all <- Rcplex(obj, constr, rhs, lb = 0, ub = 1, sense = dir,
                  objsense = "max", vtype = "B",
                  control = list(round = TRUE), n = NA)
  ## test if length of solution is correct
  ## n = 1: returns a list with 4 elements,
  ##        these are the list elements of a single solution
  ## n = 4: returns a list with 4 elements,
  ##        each element of the list represents a different solution
  ## n = NA: returns a list with >= 1 elements depending on the number
  ##        of solutions found (in this case 92)
  all(c(length(n_one) == 4, length(n_four) == 4, length(n_all) == 92) )
}
