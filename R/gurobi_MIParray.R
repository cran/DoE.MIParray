gurobi_MIParray <- function(nruns, nlevels, resolution=3, kmax=max(resolution, 2),
     distinct = TRUE, detailed = 0, start=NULL, forced=NULL, find.only=FALSE,
     maxtime = 60, nthread = 2, heurist = 0.5, MIQCPMethod=0, MIPFocus=1,
     gurobi.params=list(BestObjStop=0.5, LogFile="")){
  aufruf <- sys.call()
  nlev <- nlevels
  ## the function ensures the requested resolution and throws an error, if that is not possible
  ## the maxtime option refers to time available for for each WLP entry optimization from resolution to kmax
  ## it is in principle possible to interrupt the process and keep the result,
  ##     but unstable situations may occur, and the CPU is often not released completely

  ## check inputs
  if (!is.numeric(detailed)) stop("detailed must be numeric")
  detailed <- max(0, detailed)
  detailed <- min(3, detailed)
  detailed <- floor(detailed)  ## 0,1,2,3 possible, 1.999 -> 1 (e.g.)
  if (!is.numeric(nruns)) stop("nruns must be an integer number")
  if (!is.numeric(nlev)) stop("nlev must have integer elements")
  nfac <- length(nlev)
  if (nfac < 2) stop("nlev must have at least two elements")
  if (!is.numeric(resolution)) stop("resolution must be an integer number")
  if (!is.numeric(kmax)) stop("kmax must be an integer number")
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
  if (!length(nruns)==1) stop("nruns must be scalar")
  if (!length(resolution)==1) stop("resolution must be scalar")
  if (!length(kmax)==1) stop("kmax must be scalar")
  if (kmax < 2) stop("kmax must be at least 2")
  if (kmax > nfac) stop("kmax must not be larger than the number of factors")
  if (kmax < resolution) stop("kmax must not be smaller than resolution")
  ## added check July 2021
  if (kmax == nfac && resolution < nfac) stop("kmax must not be larger than ", nfac - 1, ", ",
               "because the last element A", nfac, " of the GWLP ",
               " is a consequence of earlier elements.")
  strength <- resolution - 1
  if (!is.logical(distinct)) stop("distinct must be logical")
  if (distinct && nruns>prod(nlev)) stop("too many runs for design with distinct runs")

  if (!is.null(start)){
    if (strength==0) stop("a start array requires resolution > 1")
    if (!is.numeric(start)) stop("start must be numeric")
    if (!all(start%%1==0)) stop("start must have integer elements")
    if (!all(start>=0)) stop("start must have non-negative elements")
    if (is.matrix(start)){
      if (!all(dim(start)==c(nruns, nfac))) stop("matrix start has wrong dimensions")
      if (!all(levels.no(start)==nlev)) stop("start and nlevels do not match")
      for (i in 1:length(nlev)) if (length(setdiff(start[,i],1:nlev[i]))>0)
        stop("invalid entries in column ", i, " of matrix start")
      start <- dToCount(start-1)
    }
    if (!length(start)==prod(nlev)) stop("vector start has wrong length")
    if (!sum(start)==nruns) stop("vector start is incompatible with nruns")
    suppressWarnings({
      if (!round(DoE.base::GWLP(countToDmixed(nlev, start), kmax=strength), 8)[strength+1]==0)
      stop("resolution of start array too low")
      })
    if (distinct && !all(start %in% c(0,1)))
      stop("start array does not comply with option distinct")
      ## fixed text in May 2021, "not" was missing
  }

  if (!is.null(forced)){
    if (strength==0) stop("forcing runs requires resolution > 1")
    if (!is.numeric(forced)) stop("forced must be numeric")
    if (!all(forced%%1==0)) stop("forced must have integer elements")
    if (!all(forced>=0)) stop("forced must have non-negative elements")
    if (is.matrix(forced)){
      if (!ncol(forced)==nfac) stop("matrix forced has wrong number of columns")
      if (nrow(forced)>=nruns) stop("matrix forced has too many rows")
      if (!all(levels.no(forced)<=nlev)) stop("forced and nlevels do not match")
      for (i in 1:length(nlev)) if (length(setdiff(forced[,i],1:nlev[i]))>0)
        stop("invalid entries in column ", i, " of matrix forced")
      forced <- dToCount(forced-1)
    }
    if (!length(forced)==prod(nlev)) stop("vector forced has wrong length")
    if (!sum(forced)<nruns) stop("vector forced fixes all runs", "use start for a start array")
    if (distinct && !all(forced %in% c(0,1)))
      stop("forced array does not comply with option distinct")
      ## fixed text in May 2021, "not" was missing
  }

  if (!is.list(gurobi.params)) stop("gurobi.params must be a named list")
  params <- gurobi.params
  if (is.null(params)) params <- vector(mode="list")
  if (is.null(params$BestObjStop)) params$BestObjStop <- 0.5
     ## BestObjStop 0.5 enforced (user can specify larger value, if desired)
     ##      for A_resolution, the bound can and will be made
     ##      using internal function lowerbounds from DoE.base
     ##      this happens in the loop
  if (params$BestObjStop < 0) {
       params$BestObjStop <- 0.5      ## objective values smaller than zero are nonsense
       warning("negative gurobi parameter BestObjStop was replaced by 0.5")
  }
  ## the strictest restriction of maximum time wins
  params$TimeLimit <- min(params$TimeLimit, maxtime)
  if (params$TimeLimit == Inf) params$TimeLimit <- NULL

  ## the largest permissible choice of Heuristics wins
  params$Heuristics <- max(params$Heuristics, heurist)
  if (params$Heuristics > 1) params$Heuristics <- 1
  if (params$Heuristics < 0) params$Heuristics <- 0

  ## the smaller thread request wins
  params$Threads <- min(nthread, params$Threads)

  if (!MIQCPMethod %in% c(-1,0,1)){
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

  ## preliminary exclusion of definitely infeasible cases
  N <- prod(nlev)
  if (strength > 0)
     if (!DoE.base::oa_feasible(nruns, nlev, strength)) stop("requested array does not exist")

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
  i <- unlist(lapply(1:N, function(obj) obj:N))
  j <- unlist(lapply(N:1, function(obj) rep(N+1-obj,obj)))
  for (kk in 1:kmax){
    Hs[[kk]] <- round(tcrossprod(mm[,(dfweg[kk]+1):dfweg[kk+1]]),8)
    Us[[kk]] <- t(mm[,(dfweg[kk]+1):dfweg[kk+1]]) ## crossprod is H
  }


  ## initialize first gurobi model (linear constraints up to strength)
  ## main effects model matrix is enforced in the linear constraint
  ##     (can be moved to being optimized by choosing resolution=1)
  qco1 <- vector(mode="list")
  ## linear constraint: sum equal to nruns
  qco1$A <- Matrix::Matrix(rep(1,N),nrow=1)              ## constrain sum
  ## orthogonal to mean all effects up to strength,
  ## cone constraint unnecessary, no additional variables

  qco1$rhs <- nruns
  qco1$sense <- "="    ## gurobi replicates it

  ## implement forced, if applicable
  if (!is.null(forced)){
    ## additional rows of A
    hilf <- which(forced>0)
    mathilf <- matrix(0, length(hilf), ncol(qco1$A))
    mathilf[,hilf] <- diag(length(hilf))  ## identity
    if (distinct) qco1 <- gurobi_modelAddLinear(qco1, mathilf, forced[hilf])
    else qco1 <- gurobi_modelAddLinear(qco1, mathilf, forced[hilf], sense=">=")
  }

  ## variable related dimension
  qco1$lb <- rep(0,N)
  qco1$ub <- rep(Inf,N)
  qco1$obj <- rep(0,N)
  if (distinct) qco1$vtype <- rep("B", N) else qco1$vtype <- rep("I", N)

  ## info elements
  qco1$info$nruns <- nruns
  qco1$info$nlev <- nlev
  qco1$info$reso <- strength+1
  qco1$info$last.k <- max(1,strength)
  qco1$info$optimizer <- "gurobi"
  qco1$info$conecur <- FALSE

  class(qco1) <- "qco"

  if (strength==0){
    qco1 <- gurobi_modelAddLinear(qco1, Us[[1]])
    qco1 <- gurobi_modelAddConeQobj(qco1, df[2])  ## for optimization of A1
  }
  else
    for (ii in 1:strength) qco1 <- gurobi_modelAddLinear(qco1, Us[[ii]])

  probs <- vector(mode="list")
  sols <- vector(mode="list")
  if (!is.null(start)){
      ## start array specified
      opt <- start
      vcur <- opt%*%Hs[[strength]]%*%opt  ## should be 0
      qco1$info$timelinear <- 0  ## no time in linear optimization
    }
  else{
      ## optimize first step
      hilf2 <- params$Heuristics
      params$Heuristics <- NULL   ## deactivate increased Heuristics for first step
      hilf3 <- params$MIPFocus
      params$MIPFocus <- NULL   ## deactivate MIPFocus for first step
      ## parameters (for documentation purpose only)
      qco1$info$params <- params
      probs[[1]] <- qco1

      qco1$info$timelinear <- system.time(sl <- gurobi::gurobi(qco1, params=params))[3]
      ## restore settings in params for further steps
      params$Heuristics <- hilf2
      params$MIPFocus <- hilf3
      qco1$info$params <- params

        sols[[1]] <- sl
      if (!sl$status == "OPTIMAL"){
        if (is.null(forced))
           stop("More runs or lower resolution needed?", " Status from initial step: ", sl$status)
        else
          stop("More runs or lower resolution needed?",
               " Or request not feasible extending from forced? ",
               " Status from initial step: ", sl$status)
      }
      opt <- sl$x     ## either A1 minimized or strength enforced
      vcur <- round(opt[1:N])%*%Hs[[max(1,strength)]]%*%round(opt[1:N])
           ## first optimum, likely 0
      ## July 2021
      stop.initial.optimal <- FALSE
      vreso <- round(opt[1:N])%*%Hs[[resolution]]%*%round(opt[1:N])  ## for checking optimality
      if (vreso <  sum(lowerbounds(nruns, nlev, resolution)) + 0.5 && kmax <= resolution){
           find.only <- TRUE
           stop.initial.optimal <- TRUE
      }
  }

  vhistory <- list(objective=vcur, counts=cbind(opt[1:N]))
  nvar <- ncol(qco1$A)

  if (strength==0) cat("=== GWLP after minimizing A1 ==========\n")
  else{
     if (is.null(start)) cat("=== GWLP after enforcing resolution ===\n")
    else cat("=== GWLP of start array ===\n")
  }
  suppressWarnings({
    print(round(DoE.base::GWLP(countToDmixed(nlev,round(opt[1:N]))),3))
  })
  cat("=======================================\n")

  ## initial step completed

  ## loop over further steps
  from <- resolution
  if (find.only) kmax <- resolution
  if (strength==0) from <- 2
  if (!kmax < from ){
    #for kmax==1 and strength==0, optimization is already finished
  for (kk in from:kmax){
    ## remove cone constraint in case of previous vcur = 0
    if (vcur==0){
      qco1 <- gurobi_modelLastQuadconToLinear(qco1)   ## should not change initial model,
                                                  ## except for strength 0
      nvar <- length(qco1$obj)
      ## increase resolution in info element, if A_R=0
      if (kk==qco1$info$reso + 1) qco1$info$reso <- kk
    }

    ## df(kk+1) +2 additional variables

    if (vcur>0){
      ## for the previous objective, the constraint to the optimum vcur is added in lb/ub
      qco1$lb[nvar-df[kk]-1] <- 0         ## t can only take integer values for integer x[1:N]
      qco1$ub[nvar-df[kk]-1] <- vcur+0.3  ## t can only take integer values for integer x[1:N]
    }

    ## now extend by current quadratic objective
    qco1 <- gurobi_modelAddLinear(qco1, Us[[kk]])
    qco1 <- gurobi_modelAddConeQobj(qco1, df[kk+1])  ## this removes any start solution

    ## the previous solution is prepared as initial solution
    ## extended by current new variables
    qco1$start <- c(opt[1:nvar],
                   round(opt[1:N])%*%Hs[[kk]]%*%round(opt[1:N]), 1,
                   Us[[kk]]%*%round(opt[1:N]))
    nvar <- as.integer(ncol(qco1$A))   #update nvar

    probs[[kk]] <- qco1
    ## July 2021
     if (find.only){
       ### not all aspects of this solution match the problem
       if (!stop.initial.optimal) sl$status <- "TIME_LIMIT" else sl$status <- "OPTIMAL"
       sl$x <- qco1$start
       sols[[kk]] <- sl
    }
    else{
    #gurobi_rsave(qco1, qco1$start)
    if (kk==qco1$info$reso) {
      ## use updated resolution, if zeroes occurred
      hilf <- params$BestObjStop
      bound <- sum(lowerbounds(nruns, nlev, kk)) + 0.5
      params$BestObjStop <- max(hilf, bound)
    }
    sl <- gurobi::gurobi(qco1, params=params)
    if (kk==qco1$info$reso){
      params$BestObjStop <- hilf
      if (!sl$status %in% c("TIME_LIMIT", "OPTIMAL") && sl$objval <= bound)
        sl$status <- "OPTIMAL"  ## better use a different word for distinguishing from Gurobi-confirmed?
    }
        sols[[kk]] <- sl
    }
    opt <- sl$x         ## optimized for A_kk
    vcur <- round(opt[1:N])%*%Hs[[kk]]%*%round(opt[1:N])  ## A_kk=vcur/nruns^2
    vhistory$objective <- c(vhistory$objective, vcur)     ## can eventually disappear
    vhistory$counts <- cbind(vhistory$counts, opt[1:N])   ## can eventually disappear
    if (!find.only){
    cat(paste("=== GWLP after optimizing A",kk," ===\n",sep=""))
      suppressWarnings({
        print(round(DoE.base::GWLP(countToDmixed(nlev,round(opt[1:N]))),3))
      })
          cat("=================================\n")
    }
  }  ## end of loop over kk
  }

  feld <- countToDmixed(nlev, round(opt[1:N])) + 1
  if (!is.null(forced)){
    ## make sure it is known which rows were forced into the array
    hilf <- countToDmixed(nlev, forced)+1
    type <- rep("new", nrow(feld))
    names(type) <- 1:nrow(feld)
    cur <- 1
    for (i in 1:nrow(feld)){
      if (all(feld[i,]==hilf[cur,])){
        cur <- cur+1
        type[i] <- "forced"
      }
      if (cur > nrow(hilf)) break
    }
  }
  class(feld) <- c("oa", "matrix")
  attr(feld, "origin") <- list(package="DoE.MIParry", call=aufruf)
  status <- lapply(sols, function(obj) obj$status)  ## list


  names(status)[1] <- paste("A1to",strength,sep="")
  ## check whether special treatment for k=1 is needed
  if (find.only)
    names(status)[resolution] <- paste0("A", resolution)
  else
    names(status)[resolution:kmax] <- paste("A",resolution:kmax,sep="")
  namdetail <- names(status)
  status <- unlist(status[sapply(status, function(obj) !is.null(obj))]) ## vector
  status[is.na(status)] <- "UNKNOWN"               ## gurobi sometimes returns NA
  laststatus <- status[length(status)]
  #vcurs <- vhistory$objective
  #names(vcurs) <- names(status)
  ausqco1 <- qco1
  ausqco1$start <- opt
  ausqco1$info$last.k <- kk
  ausqco1$info$stati <- status

  ausqco1$info$objval <- sl$objval
  ausqco1$info$objbound <- sl$objbound
  if (!is.null(forced)) ausqco1$info$forced <- type

  ## make sure info element is last
  hilf <- ausqco1$info
  ausqco1$info <- NULL
  ausqco1$info <- hilf

  if (kmax < nfac || !laststatus == "OPTIMAL")
       attr(feld, "MIPinfo") <- ausqco1
  else
       attr(feld, "MIPinfo") <- ausqco1$info ## moved here May 2021

  if (length(setdiff(status, c("OPTIMAL","USER_OBJ_LIMIT"))) > 0){
    if (laststatus=="OPTIMAL"){
      warning("Even though the last step yielded a confirmed optimum, ",
              "a previous step did not.")
    }
  }
  if (find.only) if (stop.initial.optimal) message(paste0("The initial enforcement of resolution ",
                                           "yielded the optimal A", resolution, "."))
  if (detailed) {
    names(probs) <- namdetail
    names(sols) <- namdetail
    probs[sapply(probs, is.null)] <- NULL
    sols[sapply(sols, is.null)] <- NULL
    attr(feld, "history") <- list(probs=probs, sols=sols)
  }
  if (detailed==3){
    attr(feld, "matrices") <- list(Us=Us, Hs=Hs)
  }
  feld
}

