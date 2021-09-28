## additional test cases
## these are run for avoidance of errors only
## no testoutput is created
require(DoE.MIParray)

mosek_nonverboseparams <- list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0))

if (require(gurobi)){
aus <- gurobi_MIParray(50, c(2,2,2,5,5), 2, gurobi.params = list(OutputFlag=0))
## next command must throw an error
if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
    try(mosek_MIPcontinue(aus))
    aus3improvem <- mosek_MIPcontinue(aus, improve=FALSE, maxtime=12,
                                      mosek.params = mosek_nonverboseparams)  ## fast with mosek
    aus4improvem <- mosek_MIPcontinue(aus3improvem, improve=FALSE, maxtime=12,
                                      mosek.params = mosek_nonverboseparams) ## longer time doesn't help either
  }
}
}

if (require(gurobi) && require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
    aus <- gurobi_MIParray(50, c(2,2,2,5,5), 2, find.only=TRUE, gurobi.params = list(OutputFlag=0))
    aus <- mosek_MIPcontinue(aus, maxtime=12, mosek.params = mosek_nonverboseparams)
    aus3improveg <- gurobi_MIPcontinue(aus, improve=FALSE, maxtime=20, gurobi.params = list(OutputFlag=0))  ## will not succeed to 576
    aus4improvem <- mosek_MIPcontinue(aus3improveg, improve=FALSE, maxtime=20, mosek.params = mosek_nonverboseparams)
  }
  }


## use case for which the initial array is optimal, although find.only
##   has not been specified; treated internally by setting find.only to TRUE
if (require(Rmosek))
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23")
    ausm <- mosek_MIParray(1008, c(42, 12, 6), mosek.params = mosek_nonverboseparams)
if (require(gurobi))
    ausg <- gurobi_MIParray(1008, c(42, 12, 6), gurobi.params = list(OutputFlag=0))

if (require(gurobi) && require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
    ## improve=FALSE must yield an error
    ## improve=TRUE does not make sense
    try(mosek_MIPcontinue(ausm, mosek.params = mosek_nonverboseparams))
    try(mosek_MIPcontinue(ausg, mosek.params = mosek_nonverboseparams))
    try(gurobi_MIPcontinue(ausm, gurobi.params = list(OutputFlag=0)))
    try(gurobi_MIPcontinue(ausg, gurobi.params = list(OutputFlag=0)))
  }
}

## designs found with find.only=TRUE
## make sure to leave them in the state they would be in after
## an actual optimization attempt of A_R (even if they are proven optimal)
if (require(Rmosek))
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
    ausm <- mosek_MIParray(168, c(42,12,6), 2, find.only = TRUE, mosek.params = mosek_nonverboseparams)  ## very fast
    ausmimproved <- mosek_MIPcontinue(ausm, maxtime=10, mosek.params = mosek_nonverboseparams)  ## will not improve
  }
if (require(gurobi)){
  ausg <- gurobi_MIParray(168, c(42,12,6), 2, find.only = TRUE, gurobi.params = list(OutputFlag=0))  ## very fast
  ausgimproved <- gurobi_MIPcontinue(ausg, improve=FALSE, maxtime=10, gurobi.params = list(OutputFlag=0)) ## will not improve
}
