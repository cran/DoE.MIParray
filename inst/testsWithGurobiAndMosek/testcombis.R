require(DoE.MIParray)
nruns <- 12
nlev <- c(2,2,3,2,2)

### this file tests that all combinations of situations work properly
### the settings should be such that the printed output does not depend on random inputs / computer speed etc.
###     to a reasonable extent

## mosek with resolution > 2 and kmax=resolution
if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
  aus <- mosek_MIParray(nruns,nlev, resolution=3, kmax=3,maxtime=100,
                      mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                      MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                      INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
  MIPinfo <- attr(aus, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")]) ## stable output
  }
}

if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
    ## mosek with resolution = 1 and kmax>resolution
  aus <- mosek_MIParray(nruns,nlev, resolution=1, kmax=2,maxtime=100,
                      mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                      MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                      INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
   MIPinfo <- attr(aus, "MIPinfo")
   print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output
  }
}

if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
    ## mosek with several steps of optimization also of non-zero A-values
    aus <- mosek_MIParray(nruns,nlev, resolution=2, kmax=4,maxtime=100,
                          mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                          MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                          INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
    MIPinfo <- attr(aus, "MIPinfo")
    print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output
  }}

if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
    ## mosek with resolution = 2 and kmax>resolution
    aus <- mosek_MIParray(nruns,nlev, resolution=2, kmax=3,maxtime=100,
                      mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                      MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                      INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
   MIPinfo <- attr(aus, "MIPinfo")
   print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output

   ## improve mosek solution with mosek
   ## must produce an error, because already optimal
   print(try(aus2 <- mosek_MIPcontinue(aus, improve=TRUE,
                      mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                     MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                     INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))))
## optimize next A-value with mosek
aus2 <- mosek_MIPcontinue(aus, improve=FALSE,
                          mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                          MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                          INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
MIPinfo <- attr(aus2, "MIPinfo")
print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output

  if (require(gurobi)){
    ## improve mosek solution with gurobi
    ## must produce an error, because already optimal
    print(try(aus2 <- gurobi_MIPcontinue(aus, improve=TRUE, gurobi.params = list(OutputFlag=0))))

    ## optimize next A-value with gurobi
    aus2 <- gurobi_MIPcontinue(aus, improve=FALSE, gurobi.params = list(OutputFlag=0))
    MIPinfo <- attr(aus2, "MIPinfo")
    print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output
  }
}
}

if (require(gurobi)){
  ## gurobi with resolution > 2 and kmax=resolution
  aus <- gurobi_MIParray(nruns,nlev, resolution=3, kmax=3,maxtime=100, gurobi.params = list(OutputFlag=0))
  MIPinfo <- attr(aus, "MIPinfo")
  MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")]  ## stable output
}

## gurobi with resolution = 1 and kmax>resolution
if (require(gurobi)){
  aus <- gurobi_MIParray(nruns,nlev, resolution=1, kmax=2,maxtime=100,
                      gurobi.params = list(OutputFlag=0))
  MIPinfo <- attr(aus, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output
}

## gurobi with several steps of optimization also with non-zero A-values
if (require(gurobi)){
  aus <- gurobi_MIParray(nruns,nlev, resolution=2, kmax=4,maxtime=100,
                         gurobi.params = list(OutputFlag=0))
  MIPinfo <- attr(aus, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output
}
## gurobi with resolution = 2 and kmax>resolution
if (require(gurobi)){
  aus <- gurobi_MIParray(nruns,nlev, resolution=2, kmax=3,maxtime=100,
                       gurobi.params = list(OutputFlag=0))
  MIPinfo <- attr(aus, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output

  ## improve gurobi solution with gurobi
  ## must produce an error, because already optimal
  print(try(aus2 <- gurobi_MIPcontinue(aus, improve=TRUE, gurobi.params = list(OutputFlag=0))))

  ## optimize next A-value with gurobi
  aus2 <- gurobi_MIPcontinue(aus, improve=FALSE, gurobi.params = list(OutputFlag=0))
  MIPinfo <- attr(aus2, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output

  ## improve gurobi solution with mosek
  if (require(Rmosek)){
    if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version() >= "8.1.0.23"){
      ## must produce an error, because already optimal
      print(try(aus2 <- mosek_MIPcontinue(aus, improve=TRUE,
                             mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                             MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                             INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))))

      ## optimize next A-value with mosek
      aus2 <- mosek_MIPcontinue(aus, improve=FALSE,
                                 mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                 MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                 INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))

      MIPinfo <- attr(aus2, "MIPinfo")
      print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")])  ## stable output
    }
  }
}
nruns <- 6
nlev <- c(2,2,3,2,2)
if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version() >= "8.1.0.23"){
    aus <- mosek_MIParray(nruns,nlev, resolution=2, kmax=3,maxtime=2,
                      mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                      MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                      INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
    MIPinfo <- attr(aus, "MIPinfo")
    print(MIPinfo$info[which(!names(MIPinfo$info) %in% c("timelinear", "objbound"))])  ## stable output
    ## improve mosek solution with mosek
    aus2 <- mosek_MIPcontinue(aus, mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                           MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                           INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
    MIPinfo <- attr(aus2, "MIPinfo")
    print(MIPinfo$info[which(!names(MIPinfo$info) %in% c("timelinear", "objbound"))])  ## stable output
    ## improve mosek solution with gurobi
    if (require(gurobi)){
    aus2 <- gurobi_MIPcontinue(aus, improve=TRUE, gurobi.params = list(OutputFlag=0))
    MIPinfo <- attr(aus2, "MIPinfo")
    print(MIPinfo$info[which(!names(MIPinfo$info) %in% c("timelinear", "objbound"))])  ## stable output
  }}
}

if (require(gurobi)){
    aus <- gurobi_MIParray(nruns,nlev, resolution=2, kmax=3,maxtime=2,
                           gurobi.params = list(OutputFlag=0))
    MIPinfo <- attr(aus, "MIPinfo")
    print(MIPinfo$info[which(!names(MIPinfo$info) %in% c("timelinear", "objbound"))])  ## stable output
    ## improve gurobi solution with gurobi
    aus2 <- gurobi_MIPcontinue(aus, improve=TRUE, gurobi.params = list(OutputFlag=0))
      MIPinfo <- attr(aus2, "MIPinfo")
      print(MIPinfo$info[which(!names(MIPinfo$info) %in% c("timelinear", "objbound"))])  ## stable output
    ## improve gurobi solution with mosek
      if (require(Rmosek)){
        if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version() >= "8.1.0.23"){
          aus2 <- mosek_MIPcontinue(aus, improve=TRUE,
                              mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                                                              MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                                                              INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
    MIPinfo <- attr(aus2, "MIPinfo")
    print(MIPinfo$info[which(!names(MIPinfo$info) %in% c("timelinear", "objbound"))])  ## stable output
        }
        }
}
