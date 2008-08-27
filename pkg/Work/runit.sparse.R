test.stm <- function()
  {
    stmAmat2 <- as.simple_triplet_matrix(Amat);
    checkEquals(stmAmat2,stmAmat);
  }

