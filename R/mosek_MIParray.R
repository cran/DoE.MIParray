mosek_MIParray <- function(nruns, nlevels, resolution=3, kmax=max(resolution, 2),
              distinct=TRUE, detailed=0,
              start=NULL, forced=NULL, maxtime=Inf, nthread=2,
              mosek.opts=list(verbose=10, soldetail=1),
              mosek.params=list(dparam=list(LOWER_OBJ_CUT=0.5, MIO_TOL_ABS_GAP = 0.2,
                      INTPNT_CO_TOL_PFEAS = 1e-05, INTPNT_CO_TOL_INFEAS = 1e-07),
                                iparam=list(PRESOLVE_LINDEP_USE="OFF", LOG_MIO_FREQ=100)
                                            )){
  aufruf <- sys.call()
  nlev <- nlevels
  ## INTPNT_CO_TOL_PFEAS helps, if mosek declares previous start value as infeasible
     ## default 10^-8 is too strict, in the 48 2.4.3.4.2 example even 10^-6 is too strict

  ## the function ensures the requested resolution and throws an error, if that is not possible
  ## it is in principle possible to interrupt the process and keep the result, but unstable sitautions may occur
  #kmax <- 3
  #resolution <- 3

  ## a recommendation from Romain Francois in his blog
  h <- function(w)
    if( any( grepl( "printing of extremely long output is truncated", w) ) )
    invokeRestart( "muffleWarning" )
  # withCallingHandlers( f(5), warning = h )
      ## around mosek calls will suppress annoying warnings

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
  if (!is.null(maxtime)){
    if (!is.numeric(maxtime)) stop("maxtime must be numeric")
    if (!length(maxtime)==1) stop("maxtime must be scalar")
    if (!maxtime>0) stop("maxtime must be positive")
  }
  if (!is.numeric(nthread)) stop("nthread must be numeric")
  if (!length(nthread)==1) stop("nthread must be scalar")
  if (!nthread>=0) stop("nthread must be non-negative")

  if (!length(nruns)==1) stop("nruns must be scalar")
  if (!length(resolution)==1) stop("resolution must be scalar")
  if (!length(kmax)==1) stop("kmax must be scalar")
  if (kmax < 2) stop("kmax must be at least 2")
  if (kmax > nfac) stop("kmax must not be larger than the number of factors")
  if (kmax < resolution) stop("kmax must not be smaller than the resolution")
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
      if (!all(DoE.base:::levels.no(start)==nlev)) stop("start and nlevels do not match")
      for (i in 1:length(nlev)) if (length(setdiff(start[,i],1:nlev[i]))>0)
        stop("invalid entries in column ", i, " of matrix start")
      start <- dToCount(start-1)
    }
    if (!length(start)==prod(nlev)) stop("vector start has wrong length")
    if (!sum(start)==nruns) stop("vector start is incompatible with nruns")
    if (!round(DoE.base::GWLP(countToDmixed(nlev, start), kmax=strength), 8)[strength+1]==0)
      stop("resolution of start array too low")
    if (distinct && !all(start %in% c(0,1)))
      stop("start array does comply with option distinct")
  }
  if (!is.null(forced)){
    if (strength==0) stop("forcing runs requires resolution > 1")
    if (!is.numeric(forced)) stop("forced must be numeric")
    if (!all(forced%%1==0)) stop("forced must have integer elements")
    if (!all(forced>=0)) stop("forced must have non-negative elements")
    if (is.matrix(forced)){
      if (!ncol(forced)==nfac) stop("matrix forced has wrong number of columns")
      if (nrow(forced)>=nruns) stop("matrix forced has too many rows")
      if (!all(DoE.base:::levels.no(forced)==nlev)) stop("forced and nlevels do not match")
      for (i in 1:length(nlev)) if (length(setdiff(forced[,i],1:nlev[i]))>0)
        stop("invalid entries in column ", i, " of matrix forced")
      forced <- dToCount(forced-1)
    }
    if (!length(forced)==prod(nlev)) stop("vector forced has wrong length")
    if (!sum(forced)<nruns) stop("vector forced fixes all runs", "use start for a start array")
    if (distinct && !all(forced %in% c(0,1)))
      stop("forced array does comply with option distinct")
  }

  if (!is.list(mosek.opts)) stop("mosek.opts must be a named list")
  if (!is.list(mosek.params)) stop("mosek.params must be a named list")
  if (!all(names(mosek.params) %in% c("iparam","dparam","sparam")))
      stop("invalid entries in list mosek.params")

   ## the stricter time request wins, LOWER_OBJ_CUT is always set
    ##      for A_resolution, the bound can and will be made
    ##      using hidden function lowerbounds from DoE.base
    ##      this happens in the loop
  if (maxtime<Inf){
      if (!is.null(mosek.params$dparam)){
         mosek.params$dparam$MIO_MAX_TIME <-
          min(maxtime, mosek.params$dparam$MIO_MAX_TIME)
         mosek.params$dparam$LOWER_OBJ_CUT <- max(0.5, mosek.params$dparam$LOWER_OBJ_CUT)
         }
      else mosek.params$dparam <- list(LOWER_OBJ_CUT=0.5, MIO_MAX_TIME=maxtime)
   }
   else{
     if (!is.null(mosek.params$dparam))
       mosek.params$dparam$LOWER_OBJ_CUT <- max(0.5, mosek.params$dparam$LOWER_OBJ_CUT)
     else mosek.params$dparam <- list(LOWER_OBJ_CUT=0.5)
   }

  ## the smaller thread request wins
  if (!is.null(mosek.params$iparam))
    mosek.params$iparam$NUM_THREADS <-
        min(nthread, mosek.params$iparam$NUM_THREADS)
    else mosek.params$iparam <- list(NUM_THREADS=nthread)


    opts <- mosek.opts

    ## preliminary exclusion of definitely infeasible cases
    N <- prod(nlev)
    if (strength > 0)
      if (!DoE.base::oa_feasible(nruns, nlev, strength)) stop("requested array does not exist")

  ## preparation of matrices needed in formulating the optimization problem
  df1 <- nlev-1
  D <- ff(nlev)
  df <- 1
  for (j in 1:nfac) df <- c(df,sum(apply(matrix(df1[DoE.base:::nchoosek(nfac,j)],nrow=j),2,prod)))
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
  for (kk in 1:kmax){
    Hs[[kk]] <- round(tcrossprod(mm[,(dfweg[kk]+1):dfweg[kk+1]]),8)
    Us[[kk]] <- t(mm[,(dfweg[kk]+1):dfweg[kk+1]]) ## crossprod is H
  }


    ## initialize first mosek model (linear constraints up to strength)
  ## main effects model matrix is enforced in the linear constraint
  ##     (can be moved to being optimized by choosing resolution=1)
  qco1 <- list(sense="min")
  qco1$iparam <- mosek.params$iparam
  qco1$dparam <- mosek.params$dparam
  qco1$sparam <- mosek.params$sparam
  ## linear constraint: sum equal to nruns
  qco1$A <- Matrix::Matrix(rep(1,N),nrow=1)              ## constrain sum
  qco1$bc <- matrix(nruns, 2,1)                  ## equal to nruns

  ## implement forced, if applicable
  if (!is.null(forced)){
    ## additional rows of A
    hilf <- which(forced>0)
    mathilf <- matrix(0, length(hilf), ncol(qco1$A))
    mathilf[,hilf] <- diag(length(hilf))  ## identity
    if (distinct) qco1 <- mosek_modelAddLinear(qco1, mathilf, forced[hilf])
    else qco1 <- mosek_modelAddLinear(qco1, mathilf, forced[hilf], sense=">=")
  }

  ## variable related dimension
  if (distinct)
    qco1$bx <- rbind(rep(0,N), rep(1,N))
    else qco1$bx <- rbind(rep(0,N), rep(Inf,N))         ## all nonnegative
  qco1$c <- rep(0,N)                             ## objective is constant
  qco1$intsub <- 1:N  ## give indices            ## all integer

  ## info elements
  qco1$info$nruns <- nruns
  qco1$info$nlev <- nlev
  qco1$info$last.k <- max(1,strength)
  qco1$info$optimizer <- "mosek"
  qco1$info$conecur <- FALSE
  class(qco1) <- "qco"

  ## enforced resolution by linear equality constraints
  if (strength==0){
    ## not even balance enforced, cone variables needed for A1 minimization
    qco1 <- mosek_modelAddLinear(qco1, Us[[1]])
    qco1 <- mosek_modelAddConeQobj(qco1, df[2])
  }
  else
    for (ii in 1:strength)
       qco1 <- mosek_modelAddLinear(qco1, Us[[ii]])

  probs <- vector(mode="list")
  sols <- vector(mode="list")
  if (!is.null(start)){
    ## start array specified
    opt <- start
    sl <- vector(mode="list")
    sl$sol$int$xx <- opt
    vcur <- round(opt%*%Hs[[strength]]%*%opt)  ## should be 0
    qco1$info$timelinear <- 0  ## no time in linear optimization
  }
  else{
    ## optimize first step
      probs[[1]] <- qco1
      hilf <- qco1$info
      qco1$info <- NULL
      withCallingHandlers(
        {hilf$timelinear <- system.time(sl <- Rmosek::mosek(qco1, opts=opts))[3]},
      warning = h )
      qco1$info <- hilf
      qco1$sol <- list(int=list(xx=sl$sol$int$xx))
      sols[[1]] <- sl
      if (!sl$sol$int$solsta == "INTEGER_OPTIMAL"){
           Rmosek::mosek_clean()
        if (is.null(forced))
          stop("more runs or lower resolution needed ?", " Status from initial step: problem ",
                      sl$sol$int$prosta, ", solution ", sl$sol$int$solsta)
          else
            stop("more runs or lower resolution needed ?", " Or request not feasible extending from forced? ",
                 " Status from initial step: problem ",
                 sl$sol$int$prosta, ", solution ", sl$sol$int$solsta)
           }
      opt <- sl$sol$int$xx     ## either A1 minimized or strength enforced
      vcur <- round(opt[1:N])%*%Hs[[max(1,strength)]]%*%round(opt[1:N])
           ## first optimum, likely 0
    }
  vhistory <- list(objective=vcur, counts=cbind(round(opt[1:N])))
  nvar <- ncol(qco1$A)

  if (strength==0) cat("=== GWLP after minimizing A1 ==========\n")
  else{
    if (is.null(start)) cat("=== GWLP after enforcing resolution ===\n")
    else cat("=== GWLP of start array ===\n")
  }
  print(round(DoE.base::GWLP(countToDmixed(nlev,round(opt[1:N]))),3))
  cat("=======================================\n")

  ## initial step completed

  ## loop over further steps
  from <- resolution
  if (strength==0) from <- 2
  if (from <= kmax ){
    rem <- 0  ## initialize rem
    #for kmax==1 and strength==0, optimization is already finished
  for (kk in from:kmax){
    print(kk)
    ## remove previous cone constraint in case of previous vcur = 0
    if (vcur==0){
      qco1 <- mosek_modelLastQuadconToLinear(qco1)   ## should not change initial model,
            ## except for strength 0
            ## if there was already a solution element, this is shortened
                  ## along with all x-related inputs
      ## incorporate solution vector from last optimization (which included the cone)
      rem <- length(sl$sol$int$xx)-length(qco1$c)    ## number of variables to be removed
      qco1 <- mosek_modelAddStart(qco1, sl$sol, remove=rem)
      #cat("qco1 after deleting cone variables for vcur=0\n")
      #print(lengths(qco1))
    }

    nvar <- length(qco1$c)
   if (vcur > 0){
      ## for the previous objective, the constraint to the optimum vcur is added in bx
      ## before adding the next cone
      qco1$bx[,nvar-df[kk]-1] <- c(0, vcur+0.3)
          ## t can only take integer values for integer x[1:N]
    }

    ## now extend by current quadratic objective
    ## df(kk+1) + 2 additional variables
    qco1 <- mosek_modelAddLinear(qco1, Us[[kk]])
    qco1 <- mosek_modelAddConeQobj(qco1, df[kk+1])
    if (rem > 0) sl$sol$int$xx <- sl$sol$int$xx[1:(length(sl$sol$int$xx)-rem)]
    rem <- 0  ## reset
    qco1 <- mosek_modelAddStart(qco1, sl$sol, add=Hs[[kk]])
         ## extend dimensions of sol object for new cone,
         ## including the current value of the (to be minimized) objective
    nvar <- as.integer(ncol(qco1$A))   #update nvar

    probs[[kk]] <- qco1
    #mosek_rsave(qco1, sl)
    hilf <- qco1$info
    qco1$info <- NULL
    if (kk==resolution){
      hilf2 <- qco1$dparam
      ## later, use DoE.base:::lowerbounds
      bound <- sum(lowerbounds(nruns, nlev, kk)) + 0.5
      qco1$dparam$LOWER_OBJ_CUT <- max(qco1$dparam$LOWER_OBJ_CUT, bound)
    }
    withCallingHandlers(
      {sl <- Rmosek::mosek(qco1, opts=opts)},
    warning = h )
    qco1$info <- hilf
    if (kk==resolution){
      qco1$dparam <- hilf2
      if (!sl$sol$int$solsta == "INTEGER_OPTIMAL" && sl$sol$int$pobjval <= bound)
        sl$sol$int$solsta <- "INTEGER_OPTIMAL"
          ## better use a different word for distinguishing from Mosek-confirmed?
    }
    sols[[kk]] <- sl
    opt <- sl$sol$int$xx         ## optimized for A_kk
    vcur <- round(opt[1:N])%*%Hs[[kk]]%*%round(opt[1:N])  ## A_kk=vcur/nruns^2
    vhistory$objective <- c(vhistory$objective, vcur)     ## can eventually disappear
    vhistory$counts <- cbind(vhistory$counts, opt[1:N])   ## can eventually disappear
    cat(paste("=== GWLP after optimizing A",kk," ===\n",sep=""))
    print(round(DoE.base::GWLP(countToDmixed(nlev,round(opt[1:N]))),3))
          cat("=================================\n")
  }
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
  status <- lapply(sols, function(obj) obj$sol$int$solsta)  ## list
  names(status)[1] <- paste("A1to",strength,sep="")
  namdetail <- names(status)
  names(status)[resolution:kmax] <- paste("A",resolution:kmax,sep="")
  status <- unlist(status[sapply(status, function(obj) !is.null(obj))]) ## vector
  laststatus <- status[length(status)]

  ausqco1 <- qco1
  opt[ausqco1$intsub] <- round(opt[ausqco1$intsub])
     ## necessary for the solution being accepted as starting value
  ausqco1$sol$int$xx <- opt
  ausqco1$info$last.k <- kk
  ausqco1$info$stati <- status

  ausqco1$info$objval <- sl$sol$int$pobjval
  ausqco1$info$objbound <- sl$sol$int$pobjbound
  if (!is.null(forced)) ausqco1$info$forced <- type

  if (kmax < nfac || !laststatus=="INTEGER_OPTIMAL")
       attr(feld, "MIPinfo") <- ausqco1
  if (length(setdiff(status, c("INTEGER_OPTIMAL"))) > 0){
    if (laststatus == "INTEGER_OPTIMAL"){
      warning("Even though the last step yielded a confirmed optimum, ", "a previous step did not.")
    attr(feld, "MIPinfo") <- ausqco1$info
    }
  }


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
  Rmosek::mosek_clean()  ## prevent problems that lead to machine instabilities
  feld
}

