####################################################################
## Code for handling quadratic constraints

## Constructor for building quadratic constraints
quadratic_constraint <- function( a = NULL, Q, dir = c("L", "G"), b ) {
  if( !is.null(a) )
    a <- as.double( a )
  out <- list( a = a,
               Q = as.simple_triplet_matrix(Q),
               dir = as.character(match.arg(dir)),
               b = as.double(b),
               nlinc = as.integer(length(a)),
               nqc   = as.integer(length(as.simple_triplet_matrix(Q)$v)) )
  out
  ## FIXME: currently we modify row and column indices in
  ##        simple_triplet_matrix, so that we can use them in C
  out$Q$i <- out$Q$i - 1L
  out$Q$j <- out$Q$j - 1L
  
  class(out) <- "quadratic_constraint"
  out
}

as.quadratic_constraint <- function(x, ...)
  UseMethod("as.quadratic_constraint")

as.quadratic_constraint.list <- function(x){
  quadratic_constraint(x$a, x$Q, x$dir, x$b)
}

as.quadratic_constraint.quadratic_constraint <- function(x) {
  x
}

is.quadratic_constraint <- function(x){
  inherits(x, "quadratic_constraint")
}
