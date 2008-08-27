#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <unistd.h>
#include <time.h>

#define DEFAULT_MAX_NUM_CALLS 500
#define my_error(x) {forceCplxClose = 1; error x;}

static CPXENVptr env;
static CPXLPptr lp;
static int max_numcalls;
static int numcalls = DEFAULT_MAX_NUM_CALLS;
static int forceCplxClose;

static void setparams(CPXENVptr, SEXP, int, int);
void Rcplex_init(void);
void Rcplex_close(void);
void Rcplex_free(void);
    
static SEXP getListElement(SEXP list, char *str) {
  SEXP element = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i=0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) {
      element = VECTOR_ELT(list,i);
      break;
    }
  }
  return element;
}    

void wait (int seconds) {
  clock_t endwait;
  endwait = clock() + seconds * CLOCKS_PER_SEC;
  while (clock() < endwait) {}
}

SEXP Rcplex(SEXP numcols_p,
	    SEXP numrows_p,
	    SEXP objsen_p,
	    SEXP cvec,
	    SEXP bvec,
	    SEXP Amat,
	    SEXP Qmat,
	    SEXP lb_p,
	    SEXP ub_p,
	    SEXP Rsense,
	    SEXP Rvtype,
	    SEXP isQP_p,
	    SEXP isMIP_p,
	    SEXP control)
{
  char *probname = "Rcplex";
  int numcols = INTEGER(numcols_p)[0]; 
  int numrows = INTEGER(numrows_p)[0];
  int objsen = INTEGER(objsen_p)[0];
  double *lb = REAL(lb_p); 
  double *ub = REAL(ub_p);

  char sense[numrows];
  char vtype[numcols];
  //  double dj[numrows];
  int isQP = INTEGER(isQP_p)[0];
  int isMIP = INTEGER(isMIP_p)[0];
  int trace = (int) REAL(getListElement(control,"trace"))[0];

  SEXP res;
  SEXP xopt;
  SEXP obj;
  SEXP solstat;
  SEXP extra;
  SEXP lambda;
  SEXP slack;
  SEXP nodecnt;

  int status;
  int i, j;
  int cur_numrows, cur_numcols;

  Rcplex_init();
  if(trace) Rprintf("Rcplex: num variables=%d num constraints=%d\n",numcols,numrows);

  /* lb and ub */
  for (j=0;j<numcols;++j) {
    lb[j] = R_finite(lb[j]) ? lb[j] : -CPX_INFBOUND;
    ub[j] = R_finite(ub[j]) ? ub[j] : CPX_INFBOUND;
  }

  /* set constraint inequality directions */
  for (i=0;i<numrows;++i) {
    sense[i] = CHAR(STRING_ELT(Rsense,i))[0];
  }

  /* set variable types */
  if (isMIP) {
    for (i=0;i<numcols;++i) {
      vtype[i] = CHAR(STRING_ELT(Rvtype,i))[0];
    }
  }

  /* set parameters given in control list */
  setparams(env,control,isQP,isMIP);
  lp = CPXcreateprob(env, &status, probname);

  /* check memory problems */
  if (lp == NULL) {
    my_error(("Failed to create LP.\n"));
  }

  /* copy problem data into lp */
  status = CPXcopylp(env, lp,numcols,numrows,objsen,REAL(cvec),REAL(bvec),sense,
		     INTEGER(VECTOR_ELT(Amat,0)),INTEGER(VECTOR_ELT(Amat,1)),
		     INTEGER(VECTOR_ELT(Amat,2)),REAL(VECTOR_ELT(Amat,3)),
		     lb,ub,NULL);
  if (status) {
    my_error(("Failed to copy problem data.\n"));
  }

  if (isQP) {
    status = CPXcopyquad(env,lp,
			 INTEGER(VECTOR_ELT(Qmat,0)),INTEGER(VECTOR_ELT(Qmat,1)),
			 INTEGER(VECTOR_ELT(Qmat,2)),REAL(VECTOR_ELT(Qmat,3)));
    if (status) {
      my_error(("Failed to copy quadratic term of problem data.\n"));
    }
  }

  if (isMIP) {
    status = CPXcopyctype(env,lp,vtype);
    if (status) {
      my_error(("Failed to copy vtype"));
    }
  }

  /* solve problem */
  if(isMIP) {
    status = CPXmipopt(env,lp);
  }
  else if (isQP)
    status = CPXqpopt(env,lp);
  else
    status = CPXlpopt(env,lp);

  if (status) {
    my_error(("Failed to optimize problem.\n"));
  }

  PROTECT(xopt = allocVector(REALSXP,numcols));
  PROTECT(obj = allocVector(REALSXP,1));
  PROTECT(solstat = allocVector(INTSXP,1));
  PROTECT(extra = allocVector(VECSXP,2));
  PROTECT(slack = allocVector(REALSXP,numrows));

  if (isMIP) {
    *INTEGER(solstat) = CPXgetstat(env,lp);

    status = CPXgetmipobjval(env,lp,REAL(obj));
    if (status) {
      my_error(("No MIP objective value available.\n"));
    }
    cur_numrows = CPXgetnumrows(env,lp);
    cur_numcols = CPXgetnumcols(env,lp);

    status = CPXgetmipx(env,lp,REAL(xopt),0,cur_numcols-1);
    if (status) {
      my_error(("Failed to get optimal integer x.\n"));
    }

    status = CPXgetmipslack(env,lp,REAL(slack),0,cur_numrows-1);
    if (status) {
      my_error(("Failed to get optimal slack values.\n"));
    }

    PROTECT(nodecnt = allocVector(INTSXP,1));
    *INTEGER(nodecnt) = CPXgetnodecnt(env,lp);
    SET_VECTOR_ELT(extra,0,nodecnt); SET_VECTOR_ELT(extra,1,slack);
  }
  else {
    *INTEGER(solstat) = CPXgetstat(env,lp);
    status = CPXgetobjval(env,lp,REAL(obj));
    if (status) {
      my_error(("No objective value available.\n"));
    }

    cur_numrows = CPXgetnumrows(env,lp);
    cur_numcols = CPXgetnumcols(env,lp);
    status = CPXgetx(env,lp,REAL(xopt),0,cur_numcols-1);
    if (status) {
      my_error(("Failed to get optimal x\n"));
    }

    status = CPXgetslack(env,lp,REAL(slack),0,cur_numrows-1);
    if (status) {
      my_error(("Failed to get slack values.\n"));
    }


    PROTECT(lambda = allocVector(REALSXP,numrows));
    status = CPXgetpi(env,lp,REAL(lambda),0,cur_numrows-1);
    if (status) {
      my_error(("Failed to get dual optimal dual variable values\n"));
    }
    SET_VECTOR_ELT(extra,0,lambda); SET_VECTOR_ELT(extra,1,slack);
  }
  PROTECT(res = allocVector(VECSXP,4));
  SET_VECTOR_ELT(res, 0, xopt); SET_VECTOR_ELT(res, 1, obj); SET_VECTOR_ELT(res, 2, solstat); SET_VECTOR_ELT(res,3,extra);
  UNPROTECT(7);
  status = CPXsetdefaults(env);
  if (status) {
    my_error(("Failed to set parameters to default.\n"));
  }
  return(res);
}

void Rcplex_init(void) {
  int numtries = 10;
  int status;
  char errmsg[1024];

  /* Initialize CPLEX environment */
  if (env == NULL) {
    env = CPXopenCPLEX (&status);
    while(env == NULL && numtries > 0) {
      wait(30);
      numtries--;
    }
    if (env == NULL) {
      CPXgeterrorstring(env, status, errmsg);
      error("Could not open CPLEX environment.\n%s\n",errmsg);
    }
    REprintf("CPLEX environment opened\n");
    numcalls = max_numcalls;
  } 
  else {
    numcalls--;
  }
  forceCplxClose = 0;
}

void Rcplex_close(void) {
  forceCplxClose = 1;
  Rcplex_free();
}

void Rcplex_free(void) {
  int status1, status2;
  char errmsg[1024];

  status1 = status2 = 0;
  if (lp != NULL) {
    status1 = CPXfreeprob(env,&lp);
    /*    REprintf("Freed CPLEX problem\n");*/
    lp = NULL;
  }

  if (env != NULL && (numcalls == 0 || forceCplxClose)) {
    status2 = CPXcloseCPLEX(&env);
    REprintf("Closed CPLEX environment\n");
    env = NULL;
  }

  if (status1 || status2) {
    status2 ? strcpy(errmsg,"env close ok") : CPXgeterrorstring(env,status2,errmsg);
    error("Rcplex_free failed: free problem code: %d\nClose environment msg: %s\n",status1,errmsg);
  }
  return;
}

static void setparams(CPXENVptr env, SEXP control, int isQP, int isMIP) {
  int i, status, value;
  const char *cur_parm;
  SEXP names;

  PROTECT(names = getAttrib(control,R_NamesSymbol));

  status = 1; /* avoid warning */
  for (i=0; i<length(control); i++) {
    cur_parm = CHAR(STRING_ELT(names,i));
    //    Rprintf("Cur parm: %s\n",cur_parm);

    if(strcmp(cur_parm,"trace") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_SCRIND,*REAL(VECTOR_ELT(control,i)) ? CPX_ON : CPX_OFF);
    }
    else if(strcmp(cur_parm,"method") == 0) {
      switch((value = (int) *REAL(VECTOR_ELT(control,i)))) {
      case 0:
	value = CPX_ALG_AUTOMATIC;
	break;
      case 1:
	value = CPX_ALG_PRIMAL;
	break;
      case 2:
	value = CPX_ALG_DUAL;
	break;
      case 3:
	value = CPX_ALG_NET;
	break;
      case 4:
	value = CPX_ALG_BARRIER;
	break;
      default:
	warning("Unknown optimization method %d, using default\n", value);
	value = CPX_ALG_AUTOMATIC;
      }
      if (isQP)
	status = CPXsetintparam(env,CPX_PARAM_QPMETHOD,value);
      else
	status = CPXsetintparam(env,CPX_PARAM_LPMETHOD,value);
    }
    else if(strcmp(cur_parm,"preind") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_PREIND,*REAL(VECTOR_ELT(control,i)) ? CPX_ON : CPX_OFF);
    }
    else if(strcmp(cur_parm,"aggind") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_AGGIND,(int) *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"itlim") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_ITLIM, (int) *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"epgap") == 0) {
      status = CPXsetdblparam(env,CPX_PARAM_EPGAP, *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"epagap") == 0) {
      status = CPXsetdblparam(env,CPX_PARAM_EPAGAP, *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"tilim") == 0) {
      status = CPXsetdblparam(env,CPX_PARAM_TILIM, *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"mipemphasis") == 0) {
      switch((value = (int) *REAL(VECTOR_ELT(control,i)))) {
      case 0:
	value = CPX_MIPEMPHASIS_BALANCED;
	break;
      case 1:
	value = CPX_MIPEMPHASIS_FEASIBILITY;
	break;
      case 2:
	value = CPX_MIPEMPHASIS_OPTIMALITY;
	break;
      case 3:
	value = CPX_MIPEMPHASIS_BESTBOUND;
	break;
      case 4:
	value = CPX_MIPEMPHASIS_HIDDENFEAS;
	break;
      default:
	warning("Unknown mip emphasis setting %d, using default\n", value);
	value = CPX_MIPEMPHASIS_BALANCED;
      }
      status = CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS, value);
    }
    else if(strcmp(cur_parm,"disjcuts") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_DISJCUTS,(int) *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"cliques") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_CLIQUES, (int) *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"nodesel") == 0) {
      switch((value = (int) *REAL(VECTOR_ELT(control,i)))) {
      case 0:
	value = CPX_NODESEL_DFS;
	break;
      case 1:
	value = CPX_NODESEL_BESTBOUND;
	break;
      case 2:
	value = CPX_NODESEL_BESTEST;
	break;
      case 3:
	value = CPX_NODESEL_BESTEST_ALT;
	break;
      default:
	warning("Unknown node selection strategy %d, using default\n", value);
	value = CPX_NODESEL_BESTBOUND;
      }
      status = CPXsetintparam(env,CPX_PARAM_NODESEL,value);
    }
    else if(strcmp(cur_parm,"probe") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_PROBE,(int) *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"varsel") == 0) {
      switch((value = (int) *REAL(VECTOR_ELT(control,i)))) {
      case -1:
	value = CPX_VARSEL_MININFEAS;
	break;
      case 0:
	value = CPX_VARSEL_DEFAULT;
	break;
      case 1:
	value = CPX_VARSEL_MAXINFEAS;
	break;
      case 2:
	value = CPX_VARSEL_PSEUDO;
	break;
      case 3:
	value = CPX_VARSEL_STRONG;
	break;
      case 4:
	value = CPX_VARSEL_PSEUDOREDUCED;
	break;
      default:
	warning("Unknown variable selection strategy %d, using default\n", value);
	value = CPX_VARSEL_DEFAULT;
      }
      status = CPXsetintparam(env,CPX_PARAM_VARSEL,value);
    }
    else if(strcmp(cur_parm,"flowcovers") == 0) {
      status = CPXsetintparam(env,CPX_PARAM_FLOWCOVERS,(int) *REAL(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm,"maxcalls") == 0) {
      max_numcalls = (int) *REAL(VECTOR_ELT(control,i));
    }
    else {
      warning("Unknown CPLEX parameter %s. Ignoring it.\n", cur_parm);
    }

    if (status)
      my_error(("Failure to set parameter %s, error %d.\n", cur_parm, status));
  }
  
  UNPROTECT(1);
}
