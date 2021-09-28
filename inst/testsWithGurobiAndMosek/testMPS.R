require(DoE.MIParray)
mp <- options("max.print")
options(max.print=3000)

### create reasonable scenarii
nruns <- 12
nlev <- c(2,2,3,2,2)
N <- prod(nlev)

### this file tests the functionality for exporting MPS files

  ## create list of linear problems,
  ## write out single linear problem
  ## or list of linear problems with toc
    ## list of five problems
    aus5 <- create_ILPlist(nruns,nlev)                      ## resolution 3
    aus1 <- create_ILPlist(nruns,nlev, search.orders=FALSE) ## resolution 3
    ## write a single problem
    write_MPSILPlist("test", aus1, toc=FALSE)  ## writes file test1.mps
    print(readLines("test1.MPS"), quote=FALSE)
    ## write list of problems
    write_MPSILPlist("miniproblem", aus5)
    print(probleme <- read.table("miniproblem_toc.txt"))
    print(readLines(rownames(probleme)[5]), quote=FALSE)

if (require(Rmosek)){
    if (packageVersion("Rmosek") >= "8.0.69" && Rmosek::mosek_version()>="8.1.0.23"){
  ## this part needs Rmosek because the start vector is determined by Mosek
  ## create quadratic problem with start vector
    mystart <- attr(mosek_MIParray(nruns, nlev, find.only=TRUE), "MIPinfo")$sol$int$xx[1:N]
    write_MPSMIQP("testquad", nruns, nlev, start=mystart)
  ## check also countToDmixed with new nlevels argument
    print(countToDmixed(scan("testquad.start", what=0), nlevels=nlev))
    print(readLines("testquad.MPS"), quote=FALSE)
}}
