require(DoE.MIParray)
start <- DoE.base::L16.2.8.8.1[,1:5]
force <- matrix(as.numeric(as.matrix(DoE.base::undesign(DoE.base::oa.design(L8.2.7)))), nrow=8)

## using a start value with mosek
if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
  feld <- mosek_MIParray(16, rep(2,5), resolution=4, start=start, maxtime=100,
                         mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                           MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                           INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
   print(feld)
   MIPinfo <- attr(feld, "MIPinfo")
   print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")]) ## stable output
  }}
  ## counting vector representation of the start value could also be used

dToCount(start, startfrom1 = TRUE)
  ## 32 elements for the full factorial in lexicographic order, 16 ones for the runs

## using a start value with gurobi
if (require(gurobi)){
  start <- DoE.base::L16.2.8.8.1[,1:5]
  feld <- gurobi_MIParray(16, rep(2,5), resolution=4, start=start, maxtime=100,
                         gurobi.params = list(OutputFlag=0))
  print(feld)
  MIPinfo <- attr(feld, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")]) ## stable output
}

## extending an existing array with mosek
if (require(Rmosek)){
  if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
      feld <- mosek_MIParray(16, rep(2,7), resolution=3, kmax=4, forced=force,
                    mosek.params = list(dparam=list(LOWER_OBJ_CUT = 0.5,
                        MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
                        INTPNT_CO_TOL_INFEAS = 1e-07), iparam=list(LOG=0)))
  print(feld)
  MIPinfo <- attr(feld, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")]) ## stable output
  }
}

## extending an existing array with gurobi
if (require(gurobi)){
  feld <- gurobi_MIParray(16, rep(2,7), resolution=3, kmax=4, forced=force,
                          gurobi.params = list(OutputFlag=0))
  print(feld)
  MIPinfo <- attr(feld, "MIPinfo")
  print(MIPinfo$info[which(!names(MIPinfo$info)=="timelinear")]) ## stable output
}
