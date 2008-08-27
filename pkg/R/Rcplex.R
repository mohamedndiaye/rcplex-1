# solve an lp using cplex
Rcplex <- function(cvec,Amat,bvec,Qmat=NULL,lb=0,ub=Inf,control=list(),
                   objsense=c("min","max"),sense="L",vtype=NULL)
  {
    stopifnot((is(Amat,"matrix") && is.real(Amat)) ||
              is(Amat,"dsparseMatrix") ||
              is(Amat,"simple_triplet_matrix"));
    
    numrows <- nrow(Amat);
    numcols <- ncol(Amat);

    if (!is.null(Qmat)) {
      stopifnot((is(Qmat,"matrix") && is.real(Qmat)) ||
                is(Qmat,"dsparseMatrix") ||
                is(Qmat,"simple_triplet_matrix"),
                nrow(Qmat) == numcols, ncol(Qmat) == numcols);
    }
    
    # check data dimensions
    stopifnot(length(cvec) == numcols, is.real(cvec),
              length(bvec) == numrows, is.real(bvec));


    # check bounds
    if (length(lb) == 1)
      lb <- rep(lb,numcols);

    if (length(ub) == 1)
      ub <- rep(ub,numcols);

    stopifnot(length(lb) == numcols, is.real(lb),
              length(ub) == numcols, is.real(ub))

    # check and set objective sense
    if(missing(objsense))
      objsense <- "min";
    
    stopifnot(objsense %in% c("min","max"));
    
    if (objsense == "min") {
      objsensei <- 1;
    }
    else {
      objsensei <- -1;
    }

    # check constraints sense
    if (length(sense) == 1)
      sense <- rep(sense,numrows);

    stopifnot(length(sense) == numrows, is.character(sense), all(sense %in% c('L','G','E')));

    # check variable type if needed
    if (!is.null(vtype)) {
      if (length(vtype) == 1)
        vtype <- rep(vtype,numcols);

      stopifnot(length(vtype) == numcols, is.character(vtype), all(vtype %in% c('C','I','B')));
      isMIP <- ifelse(any(vtype != "C"),1,0);
    }
    else
      isMIP <- 0;

    Acpx <- toCPXMatrix(Amat);
    Qcpx <- toCPXMatrix(Qmat);

    isQP <- ifelse(is.null(Qcpx),0,1);
    control <- check.Rcplex.control(control);
    control <- split.control.list(control);
    on.exit(.C("Rcplex_free"));
    res <- .Call("Rcplex",
                 as.integer(numcols),
                 as.integer(numrows),
                 as.integer(objsensei),
                 as.double(cvec),
                 as.double(bvec),
                 Acpx,
                 Qcpx,   
                 as.double(lb),
                 as.double(ub),
                 as.character(sense),
                 as.character(vtype),
                 as.integer(isQP),
                 as.integer(isMIP),
                 control$C);

    names(res) <- c("xopt","obj","status","extra");
    if (isMIP) {
      names(res$extra) <- c("nodecnt","slack");
      if (control$R$round) {
        intvars <- which(vtype != 'C')
        res$xopt[intvars] <- round(res$xopt[intvars])
        res$obj <- if (isQP) crossprod(res$xopt,as.matrix(Qmat)%*%res$xopt)+
          sum(cvec*res$xopt) else sum(cvec*res$xopt)
      }
    }
    else {
      names(res$extra) <- c("lambda","slack");
    }
    return(res);
  }

toCPXMatrix <- function(Amat)
  {
    if (is.null(Amat)) {
      return(NULL);
    }
    else if (is(Amat,"sparseMatrix")) {
        Amat <- as(Amat,"dgCMatrix");
        matbeg <- Amat@p;
        matcnt <- diff(c(Amat@p,length(Amat@x)));
        matind <- Amat@i;
        matval <- Amat@x;
    }
    else if (is(Amat,"simple_triplet_matrix")) {
      matbeg <- c(0L,cumsum(tabulate(Amat$j,Amat$ncol)));
      matcnt <- tabulate(Amat$j,Amat$ncol);
      matind <- Amat$i-1;
      matval <- Amat$v;
    }
    else {
      matbeg <- (0:(ncol(Amat)-1)) * nrow(Amat);
      matcnt <- rep(nrow(Amat),ncol(Amat));
      matind <- rep(0:(nrow(Amat)-1),ncol(Amat));
      matval <- as.vector(Amat);
    }
    
    return(list(matbeg=as.integer(matbeg),matcnt=as.integer(matcnt),matind=as.integer(matind),matval=as.real(matval)));
  }

check.Rcplex.control <- function(control)
  {
    con <- list(trace=1,
                maxcalls=500,
                method=0,
                preind=1,
                aggind=-1,
                itlim=1e8,
                epgap=1e-4,
                tilim=1e75,
                disjcuts=0,
                mipemphasis=0,
                cliques=0,
                nodesel=1,
                probe=0,
                varsel=0,
                flowcovers=0,
                round=0);
    
    con[names(control)] <- control;

    if (!is.null(con$trace)) {
      if(!con$trace %in% c(0,1)) {
        warning("Improper value for trace parameter. Using default.");
        con$trace <- 1;
      }
    }

    if (!is.null(con$method)) {
      if(!con$method %in% 0:4) {
        warning("Improprt value for method parameter. Using default.");
        con$method <- 0;
      }
    }
    
    if (!is.null(con$maxcalls)) {
      if(con$maxcalls <= 0) {
        warning("Improper value for maxcalls parameter. Using default.");
        con$maxcalls <- 500;
      }
    }
    if (!is.null(con$aggind)) {
      if(con$aggind < 0 && con$aggind != -1) {
        warning("Improper value for aggind parameter. Using default");
        con$aggind <- -1;
      }
    }

    if (!is.null(con$itlim)) {
      if(con$itlim < 0) {
        warning("Improper value for itlim parameter. Using default");
        con$itlim <- 1e8;
      }
    }

    if (!is.null(con$epagap)) {
      if (con$epagap) {
        warning("Improper value for epgap parameter. Using default");
        con$epagap <- 1e-6;
      }
    }
    
    if (!is.null(con$epgap)) {
      if(con$epgap < 0 || con$epgap > 1) {
        warning("Improper value for epgap parameter. Using default");
        con$epgap <- 1e-4;
      }
    }

    if (!is.null(con$tilim)) {
      if(con$tilim < 0) {
        warning("Improper value for tilim parameter. Using default");
        con$tilim <- 1e75;
      }
    }

    if (!is.null(con$mipemphasis)) {
      if(!con$mipemphasis %in% 0:4) {
        warning("Improper value for mipemphasis parameter: Using default");
        con$mipemphasis <- 0;
      }
    }

    if (!is.null(con$disjcuts)) {
      if(!con$disjcuts %in% -1:3) {
        warning("Improper value for disjcuts parameter: Using default");
        con$disjcuts <- 0;
      }
    }

    if (!is.null(con$cliques)) {
      if(!con$cliques %in% -1:2) {
        warning("Improper value for cliques parameter: Using default");
        con$cliques <- 0;
      }
    }

    if (!is.null(con$nodesel)) {
      if(!con$nodesel %in% 0:3) {
        warning("Improper value for nodesel parameter: Using default");
        con$nodesel <- 1;
      }
    }

    if (!is.null(con$probe)) {
      if(!con$probe %in% -1:3) {
        warning("Improper value for probe parameter: Using default");
        con$probe <- 0;
      }
    }

    if (!is.null(con$varsel)) {
      if(!con$varsel %in% -1:4) {
        warning("Improper value for varsel parameter: Using default");
        con$varsel <- 0;
      }
    }

    if (!is.null(con$flowcovers)) {
      if(!con$flowcovers %in% -1:2) {
        warning("Improper value for flowcovers parameter: Using default");
        con$flowcovers <- 0;
      }
    }

    if (!is.null(con$round)) {
      if (!con$round %in% c(0,1)) {
        warning("Improper value for round option: Using default");
        con$round <- 0;
      }
    }
    return(con);
  }

split.control.list <- function(control)
  {
    R.names <- c("round")
    C.names <- setdiff(names(control),R.names)
    
    R.control <- control[R.names]
    C.control <- control[C.names]
    return(list(R=R.control,C=C.control))
  }

Rcplex.close <- function() {
  invisible(.C("Rcplex_close"));
}

