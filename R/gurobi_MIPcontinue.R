gurobi_MIPcontinue <- function(qco, improve=TRUE,
     maxtime = 60, nthread = 2, heurist = 0.05, MIQCPMethod=0, MIPFocus=0,
     gurobi.params=list(BestObjStop=0.5, LogFile="")){
  aufruf <- sys.call()
  ## the function ensures the requested resolution and throws an error, if that is not possible
  ## the maxtime option refers to time available for for each WLP entry optimization from resolution to kmax
  ## it is in principle possible to interrupt the process and keep the result, but unstable sitautions may occur

  detailed <- 0 ## moved up to here July 2021
  ## check inputs
  if (!("qco" %in% class(qco) || "oa" %in% class(qco))) stop("qco must be of class qco or oa")
  if (!"qco" %in% class(qco)) {
    ## July 2021
    if ("qco" %in% class(attr(qco, "MIPinfo"))) qco <- attr(qco, "MIPinfo")
    else{
    if (!is.null(attr(qco, "history"))){
      detailed <- 1
      history <- attr(qco, "history")
      if (!is.null(attr(qco, "matrices"))){
        matrices <- attr(qco, "matrices")
        detailed <- 3
      }
    }
    ## treatment because of different find.only cases
    ## July 2021
    qco <- list(info=attr(qco, "MIPinfo"))
    }
    }
  if (is.null(qco)) stop("invalid object qco")
  if (!is.logical(improve)) stop("improve must be logical")
  nruns <- qco$info$nruns
  nlev <- qco$info$nlev
  nfac <- length(nlev)
  reso <- qco$info$reso
  last.k <- qco$info$last.k
  if (last.k==nfac-1 && !improve){
    warning("improve was set to TRUE, because A_", nfac,
            " is a consequence of earlier GWLP elements.")
    improve <- TRUE
  }
  ## refuse to improve, if already optimal
  if (improve && rev(qco$info$stati)[1] %in% c("INTEGER_OPTIMAL", "OPTIMAL"))
    stop("improvement unnecessary, ", "already optimal")
  ## append or replace optimizer info for last optimization step
  ## mosek2gurobi always replaces the last entry
  ## this needs to be changed for improve=FALSE
  if (rev(qco$info$optimizer)[1]=="mosek") {
    if (improve)
      qco <- mosek2gurobi(qco)
    else{
      hilf <- qco$info$optimizer
      if (length(hilf)==1)
        hilf <- rep(hilf, length(qco$info$stati))
      hilf <- c(hilf, "gurobi")
      qco <- mosek2gurobi(qco)
      qco$info$optimizer <- hilf
    }
  }
  if (!is.numeric(maxtime)) stop("maxtime must be numeric")
  if (!length(maxtime)==1) stop("maxtime must be scalar")
  if (!maxtime>0) stop("maxtime must be positive")
  if (!is.numeric(nthread)) stop("nthread must be numeric")
  if (!length(nthread)==1) stop("nthread must be scalar")
  if (!nthread>=0) stop("nthread must be non-negative")

  if (!is.numeric(heurist)) stop("heurist must be numeric")
  if (!length(heurist)==1) stop("heurist must be scalar")
  if (!(heurist>=0 && heurist<=1)) stop("heurist must be between 0 and 1")
  if (!is.numeric(MIQCPMethod)) stop("MIQPCMethod must be a valid integer number")

  if (!is.list(gurobi.params)) stop("gurobi.params must be a named list")
  params <- gurobi.params
  if (is.null(params)) params <- vector(mode="list")
  if (is.null(params$BestObjStop)) params$BestObjStop <- 0
     ## BestObjStop 0 enforced (user can specify larger value, if desired)
  if (params$BestObjStop < 0) {
       params$BestObjStop <- 0      ## objective values smaller than zero are nonsense
       warning("negative gurobi parameter BestObjStop was replaced by 0")
  }
  ## the strictest restriction of maximum time wins
  params$TimeLimit <- min(params$TimeLimit, maxtime)
  if (params$TimeLimit == Inf) params$TimeLimit <- NULL

  ## the larger heuristics percentage wins
  params$Heuristics <- max(params$Heuristics, heurist)
  if (params$Heuristics > 1) params$Heuristics <- 1
  if (params$Heuristics < 0) params$Heuristics <- 0

  ## the smaller thread request wins
  params$Threads <- min(nthread, params$Threads)

  if (!MIQCPMethod %in% c(0,1)){
    MIQCPMethod <- 0
    warning("invalid specification of MIQCPMethod ignored")
  }
  if (!is.null(params$MIQCPMethod)){
    if (!MIQCPMethod==params$MIQCPMethod){
      warning("conflicting specifications of MIQCPMethod ignored")
      params$MIQCPMethod <- 0  ## default
    }
  }
  else params$MIQCPMethod <- MIQCPMethod

  if (!MIPFocus %in% 0:3){
    MIPFocus <- 0
    warning("invalid specification of MIPFocus ignored")
  }
  if (!is.null(params$MIPFocus)){
    if (!MIPFocus==params$MIPFocus){
      warning("conflicting specifications of MIPFocus ignored")
      params$MIPFocus <- 0  ## default
    }
  }
  else params$MIPFocus <- MIPFocus

  N <- prod(nlev)
  kmax <- last.k + 1
  if (improve){
    ## can start immediately to try further improvement
    ## increase BestObjStop to current largest bound, if larger
    ##    (add small margin, because everything should be integer)
    ## incorporate lower bound, if last.k == reso
    if (last.k==reso) params$BestObjStop <-
        sum(lowerbounds(nruns, nlev, reso)) + 0.3
    if (qco$info$objbound > params$BestObjStop)
        params$BestObjStop <- qco$info$objbound + 0.3
    sl <- gurobi::gurobi(qco, params=params)
    opt <- sl$x  ## A_last.k improved
  }
  else{
      ## added July 2021
      if (last.k == nfac) stop("A", nfac, " is a consequence of earlier GWLP entries, ",
                                "optimizing it does not make sense.")
    ## optimize the next variant
    if (qco$info$conecur){
      nc <- length(qco$quadcon)
      poscopt <- qco$quadcon[[nc]]$Qc@i[1]+1
      start <- qco$start   ## prepare for adding this later
      vcur <- start[poscopt]  ## last objective value
      if (round(vcur,2)==0) {
         qco <- gurobi_modelLastQuadconToLinear(qco)
         if (kmax==qco$info$reso + 1) {
            qco$info$reso <- kmax
            params$BestObjStop <- max(params$BestObjStop,
                sum(lowerbounds(nruns, nlev, kmax)) + 0.5)
         }
      }
      else qco$ub[poscopt] <- vcur + 0.3   ## ensure that A_last.k cannot deteriorate
      }
    ## preparation of matrices needed in formulating the optimization problem
    df1 <- nlev-1
    D <- ff(nlev)
    df <- 1
    for (j in 1:nfac)
      df <- c(df,sum(apply(matrix(df1[nchoosek(nfac,j)],nrow=j),2,prod)))
    dfweg <- cumsum(df)
    Dfac <- as.data.frame(D)
    for (j in 1:ncol(Dfac)){
      Dfac[[j]] <- factor(Dfac[[j]])
      contrasts(Dfac[[j]]) <- DoE.base::contr.XuWu(nlev[j])
    }

    ## model matrix in normalized orthogonal coding
    mm <- model.matrix(formula(paste("~.^",kmax,sep="")), Dfac)

    ## prepare H and U matrices (Hs[[i]]==t(Us[[i]])%*%Us[[i]])
    ## Us contains transposed model matrices for main effects, 2fis, 3fis etc.
    ## Hs contains their inner products
    Us <- vector(mode="list")
    Hs <- vector(mode="list")
    ## for lower triangular representation of Hs entries
    #i <- unlist(lapply(1:N, function(obj) obj:N))
    #j <- unlist(lapply(N:1, function(obj) rep(N+1-obj,obj)))
      Hs[[kmax]] <- round(tcrossprod(mm[,(dfweg[kmax]+1):dfweg[kmax+1]]),8)
      Us[[kmax]] <- t(mm[,(dfweg[kmax]+1):dfweg[kmax+1]]) ## crossprod is H

    ## update model
        qco <- gurobi_modelAddLinear(qco, Us[[kmax]])
        qco <- gurobi_modelAddConeQobj(qco, df[kmax + 1])  ## for optimization of A_kmax
    ## add starting value
        qco$start <- c(start,
                       round(start[1:N])%*%Hs[[kmax]]%*%round(start[1:N]), 1,
                       Us[[kmax]]%*%round(start[1:N]))

    sl <- gurobi::gurobi(qco, params=params)
    opt <- sl$x     ## A_kmax minimized
    vcur <- round(opt[1:N])%*%Hs[[kmax]]%*%round(opt[1:N])
         ## next optimum
  }

    if (!improve)  cat(paste("=== GWLP after minimizing A",kmax,"===\n",sep=""))
    else            cat(paste("=== GWLP after further improving A",last.k," ==========\n",sep=""))
  suppressWarnings({
    print(round(DoE.base::GWLP(countToDmixed(nlev,round(opt[1:N]))),3))
  })
    cat("=======================================\n")

  feld <- countToDmixed(nlev, round(opt[1:N])) + 1
  class(feld) <- c("oa", "matrix")
  attr(feld, "origin") <- list(package="DoE.MIParry", call=aufruf)


  nsteps <- length(qco$info$stati)
  if (improve) qco$info$stati[nsteps] <- sl$status
  else{
    qco$info$stati <- c(qco$info$stati, sl$status)
    names(qco$info$stati)[nsteps+1] <- paste("A",kmax,sep="")
  }

  qco$start <- opt
  if (!improve) qco$info$last.k <- kmax
  if (qco$info$last.k < nfac ||
      !all(qco$info$stati %in% c("INTEGER_OPTIMAL","OPTIMAL","USER_OBJ_LIMIT")))

  qco$info$objval <- sl$objval  ## current objective value
  if (improve)
     qco$info$objbound <- max(sl$objbound, qco$info$objbound)
  else qco$info$objbound <- sl$objbound
     ## for improve, sharpest objective bound seen so far

    attr(feld, "MIPinfo") <- qco
  ## other detail
  if (detailed && improve) {
    history$sols[[length(history$sols)]] <- sl
    ## only solution changed
    ## is it beneficial to also change problem in case of changed optimizer?
    attr(feld, "history") <- history
    if (detailed==3){
      ## unchanged matrices
      attr(feld, "matrices") <- matrices
    }
  }
  if (detailed && !improve) {
    history$probs <- c(history$probs, qco)
    history$sols <- c(history$sols, sl)
    attr(feld, "history") <- history
    if (detailed==3){
      matrices$Us <- c(matrices$Us, Us[[kmax]])
      matrices$Hs <- c(matrices$Hs, Hs[[kmax]])
      attr(feld, "matrices") <- matrices
    }
  }

  feld
}

